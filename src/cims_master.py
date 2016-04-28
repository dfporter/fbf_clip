"""
CIMS master processing file.

CIMS processing pathway:
Run in the folder containing the raw fastq files in subfolders:
python ../clip/src/split_by_barcode.py

python ../clip/src/call_fastx_clipper_to_remove_adapter.py

# This requires the first 9 bases each have 20 quality, then removes them.
python ../clip/src/check_barcode_quality_and_trim.py
# The CIMS equivalent program is stripBarcode.pl,
# which can be run as:
export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/
perl ../clip/CIMS/stripBarcode.pl -len 9 -format fastq in.fastq out.fastq
# This program will filter by the barcode quality scores
# and then just call the CIMS script.

# Output files are now in ./no_barcode_plib/

python ../clip/src/call_novoalign.py

# Convert to bed format.
python ../clip/src/convert_novoalign_to_bed.py
# Calling the following:
perl ../clip/CIMS/novoalign2bed.pl --mismatch-file mm.bed in.novo out.bed

# So apparently both the mutation and tags files have the following format:
I       49342   49384   read_name  1       +
# The mutations file output by novoalign2bed.pl is in an irrelevant format.
# We have to select one type of mutation from that file and generate a
# mutations bed file with the score being the location in the read of
# that mutation type. Those files are the two arguments to CIMS.pl.
# I've written the following program to do this:

python select_mutation_dfporter.py --mut del <in.mismatch> <out.mismatch>

CIMS.pl will work correctly given this input.

"""
import glob
import sys
import os
import re
import HTSeq
import pandas
import remove_cims_with_deletions
import assign_cims_cits_to_gene
import reformat_collapsed_tags
import annotate_peaks_with_gene_type
import cims_pos_vs_fbe
from get_sequences_for_table import get_sequences_for_table

def combine_replicates(top_dir='novo_tags/'):
    top_dir = top_dir.rstrip('/') + '/'
    # Combine replicates before collapsing.
    fbf1_rep_files = []
    n2_from_fbf1_rep_files = []
    fbf2_rep_files = []
    n2_from_fbf2_rep_files = []
    for filename in glob.glob(top_dir + '/*'):
        if re.search('all', os.path.basename(filename)) is not None:
            continue
        if re.search('unknown', os.path.basename(filename)) is not None:
            continue
        if re.search('fbf', os.path.basename(filename)) is not None:
            if re.search('runx', os.path.basename(filename)) is not None:
                fbf1_rep_files.append(filename)
                print "Combining file %s with other FBF-1 replicates." % filename
            elif re.search('run813', os.path.basename(filename)) is not None:
                fbf2_rep_files.append(filename)
                print "Combining file %s with other FBF-2 replicates." % filename
        elif re.search('n2', os.path.basename(filename)) is not None:
            if re.search('runx', os.path.basename(filename)) is not None:
                n2_from_fbf1_rep_files.append(filename)
                print "Combining file %s with other N2 from FBF-1 replicates." % filename
            elif re.search('run813', os.path.basename(filename)) is not None:
                n2_from_fbf2_rep_files.append(filename)
                print "Combining file %s with other N2 from FBF-2 replicates." % filename            
    cmd = 'cat '
    for filename in fbf1_rep_files: cmd += ' %s' % filename
    cmd += '> ' + top_dir + 'all_fbf1.bed'
    print cmd
    os.system(cmd)
    cmd = 'cat '
    for filename in fbf2_rep_files: cmd += ' %s' % filename
    cmd += '> ' + top_dir + 'all_fbf2.bed'
    print cmd
    os.system(cmd)
    cmd = 'cat '
    for filename in n2_from_fbf1_rep_files: cmd += ' %s' % filename
    cmd += '> ' + top_dir + 'all_n2_from_fbf1.bed'
    print cmd
    os.system(cmd)
    cmd = 'cat '
    for filename in n2_from_fbf2_rep_files: cmd += ' %s' % filename
    cmd += '> ' + top_dir + 'all_n2_from_fbf2.bed'
    print cmd
    os.system(cmd)


def split_types():
    if not os.path.exists('mismatches_by_type'):
        os.system('mkdir mismatches_by_type')
    for filename in glob.glob('mismatch_tags/*'):
        for this_mut in ['ins', 'del', 'sub']:
            cmd = 'python ../clip/CIMS/select_mutation_dfporter.py '
            cmd += '--mut {mut} -t {tags} -i {i} -o {o}'.format(
                tags='collapsed_reformated/' + os.path.basename(filename),
                mut=this_mut, i=filename,
                o='mismatches_by_type/{base}_{mm_type}.bed'.format(
                    base=os.path.basename(filename).partition('.bed')[0],
                    mm_type=this_mut)
            )
            print cmd
            os.system(cmd)


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    #basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s


def get_sequences(bed_file, sequences,expand=0, fdr=False,
                  out_dir='cims_with_seq', highest_per_gene=False,
                  left_expand=False, right_expand=False,
                  asymmetric_expand=False,
                  ):
    print "\n*\ncims_master.get_sequences() called on %s." % bed_file
    print "\tHighest peak per gene = %s" % str(highest_per_gene)
    print "\tOutput directoy: %s" % out_dir
    if highest_per_gene: suffix = 'highest_per_gene'
    else: suffix = ''
    if not os.path.exists(out_dir):
        print "Creating output directory %s..." % out_dir
        os.system('mkdir ' + out_dir)
    if not os.path.exists('fasta'):
        os.system('mkdir fasta')
    out_filename = out_dir +'/' + os.path.basename(bed_file) + suffix
    if re.search('cims', out_dir) is None:
        fasta_filename = 'fasta/' + os.path.basename(bed_file) + 'cits' +suffix
    else:
        fasta_filename = 'fasta/' + os.path.basename(bed_file) + suffix
    fasta_file = open(fasta_filename, 'w')
    peaks = pandas.DataFrame.from_csv(bed_file, sep='\t')
    if fdr:
        peaks = peaks[peaks['FDR']<=fdr]
    if highest_per_gene:
        if 'tagNumber_k' in peaks.columns:
            remove_all_but_highest_peak_per_gene(peaks, 'tagNumber_k')
            print "\tFound k column to sort by."
        elif 'score' in peaks.columns:
            remove_all_but_highest_peak_per_gene(peaks, 'score')
            print "\tFount score column to sort by."
        else: print "Found no column to sort by."
        print "\tTried to remove all but the highest peak per gene:"
        print peaks.head()
    peaks['seq'] = ''
    for index in peaks.index:
        left = peaks.loc[index, 'left']
        right = peaks.loc[index, 'right']
        if asymmetric_expand:
            seq = sequences[peaks.loc[index, 'chrm']][
                peaks.loc[index, 'left']-left_expand:peaks.loc[index, 'right']+right_expand]
        else:
            seq = sequences[peaks.loc[index, 'chrm']][
                peaks.loc[index, 'left']-expand:peaks.loc[index, 'right']+expand]
        if peaks.loc[index, 'strand'] == '-':
            seq = rc(seq)
        peaks.loc[index, 'seq'] = seq
    peaks.to_csv(out_filename, sep='\t', index=False)
    peaks['peak_len'] = 0
    for index in peaks.index:
        peaks.loc[index, 'peak_len'] = len(peaks.loc[index, 'seq'])
    peaks = peaks[peaks['peak_len'] > 9]
    for index in peaks.index:
        if len(peaks.loc[index, 'seq']) < 40:
            continue
        fasta_file.write('>{name}\n{a_seq}\n'.format(
            name=peaks.loc[index, 'name'],
            a_seq=peaks.loc[index, 'seq']))
    fasta_file.close()


def call_cims():
    if not os.path.exists('cims_out/'):
        os.system('mkdir cims_out')
    for tag_filename in glob.glob('collapsed_reformated/*'):
        mut_filename = 'mismatches_by_type/'
        mut_filename += os.path.basename(tag_filename).partition('.bed')[0] + '_del.bed'
        out_filename = 'cims_out/' + os.path.basename(mut_filename)
        cmd = 'perl ../clip/CIMS/CIMS.pl --keep-cache -v {tags} {muts} {out}'.format(
            tags=tag_filename, muts=mut_filename, out=out_filename)
        print cmd
        os.system(cmd)




def write_fasta(peaks, filename):
    fasta_file = open(filename, 'w')
    peaks.name.apply(str)
    peaks.seq.apply(str)
    for index, row in peaks.iterrows():
        fasta_file.write('>{name}\n{a_seq}\n'.format(
            name=peaks.loc[index, 'name'],
            a_seq=peaks.loc[index, 'seq']))
    fasta_file.close()

def get_sequences_for_all_cims_out(
    expand=20, top_dir='cims_out', highest_per_gene=False,
    out_dir='cims_with_seq/', fdr=0.001):
    """Add a sequence column, and write a fasta file.
    """
    fasta_filename = '/scratch/indexes/WS235.fa'
    sequences = dict(
        (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for cims_file in glob.glob(top_dir + '/*'):
        get_sequences(cims_file, sequences, expand=expand, fdr=fdr,
                      highest_per_gene=highest_per_gene,
                      out_dir=out_dir)


def cits():
    # CITS analysis.
    # 1. Remove reads with CIMS.
    # 2. Run bedExt.pl to shift the bed file coordinates (both start and end)
    # one nt to the left.
    # 3. awk
    # 4. Get peaks with tag2peak.pl. This is the CITS result.
    if not os.path.exists('shifted_novo_tags/'):
        os.system('mkdir shifted_novo_tags')
    if not os.path.exists('shifted_clusters/'):
        os.system('mkdir shifted_clusters')
    if not os.path.exists('shifted_clean/'):
        os.system('mkdir shifted_clean')
    if not os.path.exists('cits_out/'):
        os.system('mkdir cits_out')
    if not os.path.exists('cache/'):
        os.system('mkdir cache')
    for bed_filename in glob.glob('collapsed_reformated/*.bed'):  #runx_fbf_tggc_adapter_removed_no_barcode.bed'):
        print """
*\n*\n*\n\n*\n*\n\n*\n*\n\n*\n*\n\n*\n*\n\n*\n*\n\n*\n*\n
"""
        os.system('rm cache/*')
        shifted_bed_filename = 'shifted_novo_tags/' + os.path.basename(bed_filename)
        if os.path.exists(shifted_bed_filename):
            # bedExt.pl will not overwrite existing files.
            os.system('rm ' + shifted_bed_filename)
        cmd = 'perl ../clip/CIMS/bedExt.pl -n up -l "-1" -r "-1" -v'
        cmd += ' {infile} {outfile}'.format(
            infile=bed_filename, outfile=shifted_bed_filename)
        print cmd
        os.system(cmd)
        # With the current method of collapsing reads, no names will match
        # with the mismatches_by_type file because [][][] exist in name.
        remove_cims_with_deletions.remove_deletions(
            shifted_filename=shifted_bed_filename)
        # Output is in in ./no_deletions/
        # Cluster.
        clustered_bed_filename = 'shifted_clusters/' + os.path.basename(
            bed_filename)
        clustered_bed_filename_zero = clustered_bed_filename + '.0'
        # Unclear how this could work with maxgap -1, which is
        # protocol, since bedExt.pl reduces you to a single nucleotide.
        # I get no clusters with maxgap -1, and do get clusters with 10.
        cmd = 'perl ../clip/CIMS/tag2cluster.pl -s -maxgap "-1" -of bed -v'
        cmd += ' {infile} {outfile}'.format(
            infile=bed_filename,  # Should this be no_deletions/base?
            outfile=clustered_bed_filename_zero)
        print cmd
        os.system(cmd)
        cmd = """awk '{if($5>2) {print $0}}' %s > %s""" % (
            clustered_bed_filename_zero, clustered_bed_filename)
        print cmd
        os.system(cmd)
        # Fix clean_filename to be reads with no deletions.
        clean_filename= 'no_deletions/' + os.path.basename(bed_filename)
        cits_out_filename = 'cits_out/' + os.path.basename(bed_filename)
        cmd = 'perl ../clip/CIMS/tag2peak.pl --keep-cache -c ./cache -ss -v -gap 25'
        cmd += ' -p 0.001 {cluster} {no_cims} {sthirty}'.format(
            cluster=clustered_bed_filename, no_cims=clean_filename,
            sthirty=cits_out_filename)
        print cmd
        os.system(cmd)


def remove_all_but_highest_peak_per_gene(peaks, colname):
    #gene_ranks = append(peaks.loc[:, ('gene_name', 'height')])
    peaks.sort(colname, inplace=True, ascending=False)
    peaks.drop_duplicates('gene_name', inplace=True)


def apply_filter(peaks, cits=False, fdr=False,
                 min_k=False, min_score=False,
                 highest_per_gene=False):
    """Applies a filter to dataframe.
    """
    if fdr and not cits:
        peaks = peaks[peaks['FDR'] < fdr]
    if min_k and not cits:
        peaks = peaks[peaks['tagNumber_k']>=min_k]
    return peaks


def subset_peaks_in_folder_that_are_reproducible(input_folder,
                                                 label='cims'):
#    for filename in glob.glob(input_folder + '/*'):
    subset_peaks_in_a_file_by_rep(input_folder, 'fbf1',
                                  min_reps=2, label=label)
    subset_peaks_in_a_file_by_rep(input_folder, 'fbf2',
                                  min_reps=2, label=label)
    subset_peaks_in_a_file_by_rep(input_folder, 'both',
                                  min_reps=5, label=label)
    


def subset_peaks_in_a_file_by_rep(input_folder, fbf, min_reps=3, label='cims'):
    if fbf == 'both':
        rep1 = pandas.read_csv(input_folder + '/fbf1_gcca.txt', sep='\t')
        rep2 = pandas.read_csv(input_folder + '/fbf1_tccg.txt', sep='\t')
        rep3 = pandas.read_csv(input_folder + '/fbf1_aacc.txt', sep='\t')
        rep4 = pandas.read_csv(input_folder + '/fbf2_gcca.txt', sep='\t')
        rep5 = pandas.read_csv(input_folder + '/fbf2_tccg.txt', sep='\t')
        rep6 = pandas.read_csv(input_folder + '/fbf2_aacc.txt', sep='\t')
        reps_df_list = [rep1, rep2, rep3, rep4, rep5, rep6]
    else:
        rep1 = pandas.read_csv(input_folder + '/' + fbf + '_gcca.txt', sep='\t')
        rep2 = pandas.read_csv(input_folder + '/' + fbf + '_tccg.txt', sep='\t')
        rep3 = pandas.read_csv(input_folder + '/' + fbf + '_aacc.txt', sep='\t')
        reps_df_list = [rep1, rep2, rep3]
#    all_reps = pandas.read_csv(
#        input_folder + '/combined_' + fbf + '.txt', sep='\t')
    cims_sets = load_replicates(reps_df_list)
    all_pos_count = {}
    for rep_num in cims_sets:
        for pos in cims_sets[rep_num]:
            all_pos_count.setdefault(pos, 0)
            all_pos_count[pos] += 1
    for rep in reps_df_list:
        rep['in_n_reps'] = 0
        rep['loc_str'] = ''
        rep.name.apply(str)
        rep.chrm.apply(str)
        rep.left.apply(int)
        for index, row in rep.iterrows():
            tup = (rep.loc[index, 'chrm'], rep.loc[index, 'left'])
            rep.loc[index, 'in_n_reps'] = all_pos_count[tup]
            rep.loc[index, 'loc_str'] = rep.loc[index, 'chrm'] + '_' + str(rep.loc[index, 'left'])
    all_reps = pandas.DataFrame(columns=reps_df_list[0].columns)
    num_repod = 0
    for rep in reps_df_list:
        for index, row in rep.iterrows():
            if rep.loc[index, 'in_n_reps'] < min_reps:
                continue
            if rep.loc[index, 'loc_str'] in all_reps['loc_str']:
                i = all_reps[all_reps['loc_str']==rep.loc[index, 'loc_str']].index[0]
                if rep.loc[index, 'height'] > all_reps.loc[i, 'height']:
                    all_reps.loc[i] = rep.loc[index]
            else:
                all_reps.loc[num_repod] = rep.loc[index]
                num_repod += 1
    all_reps.name.apply(str)
    all_reps.seq.apply(str)
#    for index, row in all_reps.iterrows():
#        all_reps.loc[index, 'loc_str'] = all_reps.loc[index, 'chrm'] + '_' + str(all_reps.loc[index, 'left'])
    all_reps.sort(['score'], ascending=False, inplace=True)
#    all_reps.drop_duplicates('loc_str', inplace=True)
    os.system('mkdir reproducible')
    all_reps.to_csv(
        'reproducible/combined_' + fbf + '_' + label + '.txt', sep='\t',
        index=False)
    write_fasta(all_reps, 'fasta/reproducible' + fbf + '_' + label + '.txt')


def load_replicates(list_of_dataframes):
    cims_sets = {}
    for num, a_dataframe in enumerate(list_of_dataframes):
        cims_sets[num] = set()
        for index in a_dataframe.index:
            cims_sets[num].add(
                (a_dataframe.loc[index, 'chrm'],
                 a_dataframe.loc[index, 'left']))
    return cims_sets
#peaks = pandas.read_csv('reproducible/combined_fbf1.txt', sep='\t')
#fasta_filename = '/scratch/indexes/WS235.fa'
#sequences = dict(
#    (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
#write_fasta(peaks, 'fasta/reproducible_fbf1.fa')

ipython = '''
import pandas
peaks = pandas.read_csv('cims_tables_one_per_gene/combined_fbf1.txt', sep='\t')
peaks['n_fbe'] = 0
for index, row in peaks.iterrows():
    peaks.loc[index, 'seq'] = peaks.loc[index, 'seq'].lower()
    peaks.loc[index, 'n_fbe'] = len(
        re.findall('tgt\w\w\wat', peaks.loc[index, 'seq']))
'''

def adjust_peak_width(input_folder, table_dir='cims_alt_tables/'):
    fasta_filename = '/scratch/indexes/WS235.fa'
    if not os.path.exists(table_dir):
        os.system('mkdir ' + table_dir)
    sequences = dict(
        (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for filename in glob.glob(input_folder + '/*'):
        peaks = pandas.read_csv(filename, sep='\t')
        print peaks.head()
        get_sequences_for_table(peaks, sequences,
                                expand=10, max_width=False)
        write_fasta(peaks, 'fasta/' + os.path.basename(filename))
        peaks.sort('height', ascending=0, inplace=True)
        write_peaks_table(peaks, filename, tables_folder=table_dir)


def rename_columns(peaks, cits=False):
    print peaks.head()
    if cits:
        peaks.columns = ['chrm', 'left', 'right', 'name',
                         'score', 'strand']
        peaks['height'] = peaks['score']
    else:
        peaks.columns = ['chrm', 'left', 'right', 'name', 'score', 'strand',
                         'tagNumber_k', 'mutationFreq_m', 'FDR',
                         'count_at_least_m_k']
        peaks['height'] = peaks['tagNumber_k']

def write_peaks_table(peaks, filename, tables_folder='cims_tables/'):
    if not os.path.exists(tables_folder):
        os.system('mkdir ' + tables_folder)
    if re.search('runx', os.path.basename(filename)) is not None:
        out_dir = 'fbf1/'
    elif re.search('fbf1', os.path.basename(filename)) is not None:
        out_dir = 'fbf1/'
    else: out_dir = 'fbf2/'
    out_name = False
    if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
    if (re.search('fbf_cgga', os.path.basename(filename)) is not None) or (
        re.search('fbf_barcode3', os.path.basename(filename)) is not None):
        if out_dir == 'fbf1/':
            out_name = 'fbf1_tccg.txt'
        else: out_name = 'fbf2_tccg.txt'
    if (re.search('fbf_ggtt', os.path.basename(filename)) is not None) or (
        re.search('fbf_barcode2', os.path.basename(filename)) is not None):
        if out_dir == 'fbf1/':
            out_name = 'fbf1_aacc.txt'
        else: out_name = 'fbf2_aacc.txt'
    if (re.search('fbf_tggc', os.path.basename(filename)) is not None) or (
        re.search('fbf_barcode1', os.path.basename(filename)) is not None):
        if out_dir == 'fbf1/':
            out_name = 'fbf1_gcca.txt'
        else: out_name = 'fbf2_gcca.txt'
    if re.search('all', os.path.basename(filename)) is not None:
        if out_dir == 'fbf1/':
            out_name = 'combined_fbf1.txt'
        else: out_name = 'combined_fbf2.txt'
    if out_name:
        peaks.to_csv(tables_folder + '/%s' % out_name,
                     sep='\t', index=False)
    else:
        peaks.to_csv(tables_folder + '/' + os.path.basename(filename),
                     sep='\t', index=False)


def create_tables(cims=True):
    # Load some library files.
    print "create_tables() called."
    fasta_filename = '/scratch/indexes/WS235.fa'
    sequences = dict(
        (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
    gtf_df = pandas.read_csv('../clip/lib/gtf_with_biotype_column.txt', sep='\t')
    if cims:
        globstr = 'cims_out/*'
        cits_option = False
        fdr = 0.001
        table_dir = 'cims_tables/'
    else:
        create_cits_tables(sequences, gtf_df)
        return
    for filename in glob.glob(globstr):
        # Apply filter also renames columns.
        print "create_tables(): %s" % filename
        print "Loading peaks..."
        peaks = pandas.read_csv(filename, sep='\t')
        #peaks = peaks.head()
        rename_columns(peaks)
        peaks = apply_filter(peaks, cits=cits_option, fdr=fdr)
        assign_cims_cits_to_gene.assign_table(
            peaks, gtf_df=gtf_df, given_gtf=True)
        annotate_peaks_with_gene_type.annotate_peaks_with_gene_type(
            peaks, gtf_filename='../clip/lib/gtf_with_biotype_column.txt')
        peaks = peaks[peaks['biotype']=='protein_coding']
        get_sequences_for_table(peaks, sequences, expand=10)
        write_fasta(peaks, 'fasta/' + os.path.basename(filename))
        peaks.sort('height', ascending=0, inplace=True)
        write_peaks_table(peaks, filename, tables_folder=table_dir)
        print "Finished processing..."


def create_cits_tables(sequences, gtf_df):
    globstr = 'cits_out/*'
    cits_option = True
    p = 0.001
    table_dir = 'cits_tables/'
    if not os.path.exists(table_dir):
        os.system('mkdir ' + table_dir)
    for filename in glob.glob(globstr):
        print "create_tables(): %s" % filename
        peaks = pandas.read_csv(filename, sep='\t', header=None)
        rename_columns(peaks, cits=True)
        assign_cims_cits_to_gene.assign_table(
            peaks, gtf_df=gtf_df, given_gtf=True)
        annotate_peaks_with_gene_type.annotate_peaks_with_gene_type(
            peaks, gtf_filename='../clip/lib/gtf_with_biotype_column.txt')
        peaks = peaks[peaks['biotype']=='protein_coding']
        get_sequences_for_table(peaks, sequences, expand=20,
                                asymmetric_expand=False,
                                left_expand=20,
                                right_expand=20)
        write_fasta(peaks, 'fasta/' + os.path.basename(filename).partition(
            '.bed')[0] + '_cits.fa')
        write_peaks_table(peaks, os.path.basename(filename),
                          tables_folder=table_dir)


def create_folder_of_highest_per_gene(
            in_folder='cims_tables/',
            out_folder='cims_tables_one_per_gene/',
            cims=True):
    if not os.path.exists(in_folder):
        os.system('mkdir ' + in_folder)
    if not os.path.exists(out_folder):
        os.system('mkdir ' + out_folder)
    for filename in glob.glob(in_folder + '/*'):
        print "Outputing the highest peak per gene for %s..." % filename
        peaks = pandas.read_csv(filename, sep='\t')
        remove_all_but_highest_peak_per_gene(peaks, 'height')
        peaks.to_csv(out_folder + '/' + str(
            os.path.basename(filename)), sep='\t', index=False)
        write_fasta(peaks, 'fasta/' + os.path.basename(filename) \
                    + 'highest_per_gene.fa')


def score_positives(peaks):
    known_pos = set(['gld-1', 'htp-1', 'htp-2', 'mpk-1', 'him-3',
                         'fbf-1', 'lip-1', 'syp-2', 'fbf-2', 'fog-1',
                         'fem-3', 'syp-3', 'gld-3', 'fog-3', 
                         'egl-4'])
    try:
        obs_genes = set(peaks['gene_name'])
    except:
        print "No gene_name column in peaks table: %s" % str(peaks.head())
        return {'observed positives': set(),
                'number of observed positives': 0,
            'missing positives': list(known_pos),
                'number of missing positives': len(list(known_pos)),
            'expected': len(list(known_pos))
                }
    obs_pos = known_pos & obs_genes
    missing_pos = known_pos - obs_genes
    obs_pos_n = len(list(obs_pos))
    missing_pos_n = len(list(missing_pos))
    return {'observed positives': obs_pos, 'number of observed positives': obs_pos_n,
            'missing positives': missing_pos, 'number of missing positives': missing_pos_n,
            'expected': len(list(known_pos))}


def score_metric(filename, label=""):
    if not label:
        label = os.path.basename(filename)
    peaks = pandas.read_csv(filename, sep='\t')
    if 'chrm' not in peaks.columns:
        print "Incorrect or absent header - no chrm."
        return "NA"
    if len(peaks['chrm']) == 0:
        print "score_binding_site given a dataframe without peaks."
        return "NA"
    score_binding_site(peaks)
    positives = score_positives(peaks)
    return write_metrics(peaks, positives, label)


def score_binding_site(peaks):
    peaks.seq.apply(str)
    peaks['has_fbe'] = 0
    if len(peaks['chrm']) == 0:
        print "score_binding_site given a dataframe without peaks."
    for index, row in peaks.iterrows():
        if not isinstance(peaks.loc[index, 'seq'], basestring):
            peaks.loc[index, 'has_fbe'] = 0
            continue
        if re.search('tgt\w\w\wat', peaks.loc[index, 'seq'].lower()) is not None:
            peaks.loc[index, 'has_fbe'] = 1
        else:
            peaks.loc[index, 'has_fbe'] = 0


def write_metrics(peaks, positives, label):
    li = """
Dataset: {label}
Number of peaks: {df_size}
Number of genes: {n_genes}
With FBE: {with_fbe}, {fbe_perc}%
Without FBE: {without_fbe}
Positive controls: {observed}/{expected}
Missing positives: {missing}
""".format(label=label,
    df_size=len(peaks), n_genes=len(list(set(peaks['gene_name']))),
           with_fbe=len(peaks[peaks['has_fbe']==1]),
           fbe_perc= float(100 * len(peaks[peaks['has_fbe']==1])/len(peaks)),
           without_fbe=len(peaks[peaks['has_fbe']==0]),
           observed=positives['number of observed positives'],
           expected=positives['expected'],
           missing=positives['missing positives'])
    print li
    return li

if __name__ == '__main__':
    for_cut_and_paste="""
export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/
    """
    cmd = r'''export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/'''
    print cmd
    #os.system(cmd)

    cmd = 'python ../clip/src/split_by_barcode.py'
    #print cmd
    #os.system(cmd)

    cmd = 'python ../clip/src/call_fastx_clipper_to_remove_adapter.py fbf-1andN2/'
    #print cmd
    #os.system(cmd)
    cmd = 'python ../clip/src/call_fastx_clipper_to_remove_adapter.py fbf-2andN2/'
    #print cmd
    #os.system(cmd)

    cmd = 'python ../clip/src/check_barcode_quality_and_trim.py'
    print cmd
    #os.system(cmd)

    cmd = 'python ../clip/src/call_novoalign.py'
    print cmd
    #os.system(cmd)

    cmd = 'python ../clip/src/convert_novoalign_to_bed.py'
    print cmd
    #os.system(cmd)

    print "Combining replicates in novo_tags/..."
    #combine_replicates(top_dir='novo_tags/')

    cmd = 'python ../clip/src/collapse_duplicate_tags_for_cims.py'
    print cmd
    #os.system(cmd)

#    reformat_collapsed_tags.reformat_collapsed_tags(
#        input_dir='novo_tags_collapse/',
#        output_dir='collapsed_reformated/')
#    split_types()
    #sys.exit()
    #call_cims()  # Writes to cims_out/
    #cits()  # Writes to cits_out/
#    create_tables(cims=False)  # Writes to cits_tables/
#    sys.exit()
#    create_tables(cims=True)  # Writes to cims_tables/
#    create_folder_of_highest_per_gene(
#        in_folder='cims_tables/',
#        out_folder='cims_tables_one_per_gene/',
#        cims=True)
#    adjust_peak_width('cims_tables/', table_dir='cims_10_width/')
#    subset_peaks_in_folder_that_are_reproducible('cims_tables/')
#    subset_peaks_in_folder_that_are_reproducible('cits_tables/',
#                                                 label='cits')
#    subset_peaks_in_folder_that_are_reproducible('cims_tables/',
#                                                 label='cims')
#    for filename in glob.glob('cims_tables/*'):
#        print "\n*****\ncims_tables:"
#        score_metric(filename, label='')
#    for filename in glob.glob('cits_tables/*'):
#        print "\n*****\ncits_tables:"
#        score_metric(filename, label='')
#    for filename in glob.glob('reproducible/*'):
#        print "\n*****\nreproducible:"
#        score_metric(filename, label='')
    cims_pos_vs_fbe.cim_vs_site_in_dir('reproducible/')
