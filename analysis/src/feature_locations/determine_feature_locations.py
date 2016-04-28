"""
Determines the average locations of an FBF peak, FBE, poly(A),
txpt bounds and CDS bounds.

1. Loads the entire gtf.
2. Load the peaks.
3. Subset the gtf to the targets.

Feature locations will be held in a class, flocs.
4. Marks the CDS and txpt bounds with mark_txpts.
"""
import HTSeq
from flocs import flocs
from build_image_objects_for_heatmaps import *
import peak_in_gene_line_plot
import scatterplot_correlation_by_wig
import output_heatmap_figures
import logging
import datetime

norm_left_len = 200.
norm_right_len = 200.
norm_cds_len = 1000.
utr_len = 800
verbose = False

logger = logging.getLogger(__name__)

def load_gtf(in_filename, return_split_lines=False):
    utrs = {}
    exon = {}
    cds = {}
    other = {}
    print "load_gtf(): Function called on %s as input." % in_filename
    lt = pandas.read_csv(in_filename, sep='\t')
    with open(in_filename, 'r') as f:
        header = next(f)
        #print "saved header %s" % header
        for li in f:
            #if not re.search('ced-12', li): continue
            s = li.rstrip('\n').split('\t')
            # biotype = s[13]
            #if s[13] != 'protein_coding': continue 
            #name = s[9]
            #txpt = s[10]
            if s[2] == 'UTR':
                utrs.setdefault(s[9], {})  # Gene name first.
                utrs[s[9]].setdefault(s[10], [])  # Then txpt id.
                if return_split_lines:
                    utrs[s[9]][s[10]].append(li.rstrip('\n').split('\t'))
                    if len(utrs[s[9]][s[10]][-1]) < 5:
                        print "wrong len: %s" % str(utrs[s[9]][s[10]][-1])
                    utrs[s[9]][s[10]][-1][3] = int(utrs[s[9]][s[10]][-1][3])
                    utrs[s[9]][s[10]][-1][4] = int(utrs[s[9]][s[10]][-1][4])
                else:
                    utrs[s[9]][s[10]].append(li)
            elif s[2] == 'exon':
                exon.setdefault(s[9], {})  # Gene name first.
                exon[s[9]].setdefault(s[10], [])  # Then txpt id.
                if return_split_lines:
                    exon[s[9]][s[10]].append(li.rstrip('\n').split('\t'))
                    exon[s[9]][s[10]][-1][3] = int(exon[s[9]][s[10]][-1][3])
                    exon[s[9]][s[10]][-1][4] = int(exon[s[9]][s[10]][-1][4])
                else:
                    exon[s[9]][s[10]].append(li)
            elif s[2] == 'CDS':
                cds.setdefault(s[9], {})  # Gene name first.
                cds[s[9]].setdefault(s[10], [])  # Then txpt id.
                if return_split_lines:
                    cds[s[9]][s[10]].append(li.rstrip('\n').split('\t'))
                    cds[s[9]][s[10]][-1][3] = int(cds[s[9]][s[10]][-1][3])
                    cds[s[9]][s[10]][-1][4] = int(cds[s[9]][s[10]][-1][4])
                else:
                    cds[s[9]][s[10]].append(li)
            else:
                other.setdefault(s[9], {})  # Gene name first.
                other[s[9]].setdefault(s[10], [])  # Then txpt id.
                if return_split_lines:
                    other[s[9]][s[10]].append(li.rstrip('\n').split('\t'))
                    other[s[9]][s[10]][-1][3] = int(other[s[9]][s[10]][-1][3])
                    other[s[9]][s[10]][-1][4] = int(other[s[9]][s[10]][-1][4])
                else:
                    other[s[9]][s[10]].append(li)
    print "load_gtf(): created dicts of genes size: utrs: %i exons: %i cds: %i other: %i" % (
        len(utrs), len(exon), len(cds), len(other))
    print "load_gtf(): total txpt numbers are %i %i %i %i" %(
        sum([len(utrs[gene_name] )for gene_name in utrs]),
        sum([len(exon[gene_name]) for gene_name in exon]),
        sum([len(cds[gene_name]) for gene_name in cds]),
        sum([len(other[gene_name]) for gene_name in other])
        )
    return (utrs, exon, cds, other)


def create_file_of_longest_transcript_per_gene(
        in_filename='', out_filename='../clip/lib/gtf_one_txpt_per_gene.txt'):
    (utrs, exon, cds, other) = load_gtf(in_filename)
    longest_txpt_per_gene = {}
    without_utr = 0
    for gene in cds:
        txpt_len = []
        if (gene not in utrs) or (len(utrs[gene]) < 2): 
            without_utr += 1
            for txpt in exon[gene]:  # All CDS genes have an exon line.
                _left = [(int(ex.split('\t')[3]), txpt) for ex in exon[gene][txpt]]
                _right = [(int(ex.split('\t')[4]), txpt) for ex in exon[gene][txpt]]
                _left = sorted(_left, key=lambda x: x[0])
                _right = sorted(_right, key=lambda x: x[0])
                txpt_len.append((_right[-1][0] - _left[0][0], txpt))               
        else:
            for txpt in utrs[gene]:
                left_utrs = [(int(a_utr.split('\t')[3]), txpt) for a_utr in utrs[gene][txpt]]
                right_utrs = [(int(a_utr.split('\t')[4]), txpt) for a_utr in utrs[gene][txpt]]
                left_utrs = sorted(left_utrs, key=lambda x: x[0])
                right_utrs = sorted(right_utrs, key=lambda x: x[0])
                txpt_len.append((right_utrs[-1][0] - left_utrs[0][0], txpt))
        longest_txpt_per_gene[gene] = sorted(txpt_len, key=lambda x: x[0])[-1][1]
    print "create_file_of_longest_transcript: mapped %i genes to transcripts." % len(longest_txpt_per_gene)
    print "create_file_of_longest_transcript: example: %s->%s" % (
        longest_txpt_per_gene.keys()[0], longest_txpt_per_gene.values()[0])
    print "create_file_of_longest_transcript: genes with exons but < 2 utrs: %i" % (
        without_utr)
    longest_txpt_set = set(longest_txpt_per_gene.values())
    with open(in_filename, 'r') as f:
        header = next(f)
    utr_genes = set(utrs.keys())
    cds_genes = set(cds.keys())
    exon_genes = set(exon.keys())
    cds_no_exon = cds_genes - exon_genes
    exon_no_cds = exon_genes - cds_genes
    exon_and_cds = exon_genes & cds_genes
    print """create_file_of_longest_transcript: \
cds, no exon: %i
exon, no cds: %i
exon and cds: %i""" % (
    len(list(cds_no_exon)), len(list(exon_no_cds)), len(list(exon_and_cds)))
    with open(out_filename, 'w') as f:
        f.write(header)
        for gene in utrs:
            for txpt in utrs[gene]:
                if txpt in longest_txpt_set:
                    for a_utr in utrs[gene][txpt]:
                        f.write(a_utr)
        for gene in exon:
            for txpt in exon[gene]:
                if txpt in longest_txpt_set:
                    for a_exon in exon[gene][txpt]:
                        f.write(a_exon)
        for gene in cds:
            for txpt in cds[gene]:
                if txpt in longest_txpt_set:
                    for a_cds in cds[gene][txpt]:
                        f.write(a_cds)
        for gene in other:
            for txpt in other[gene]:
                if txpt in longest_txpt_set:
                    for a_exon in other[gene][txpt]:
                        f.write(a_exon)


def mark_txpts(sequences, rc_sequences, peaks, utrs, exon, cds, other):
    if verbose: print "mark_txpts(): Adding gene sequences..."
    # This should be txpt_id. Fix later.
#    for txpt_id in set(gtf_df['transcript_id'].tolist()):
    flocs_map = {}
    for n, gene_name in enumerate(set(peaks['gene_name'].tolist())):
        if verbose:
            print str("*" * 14) + "\nmark_txpts(): Marking features in gene_name %s" % gene_name
        if gene_name not in cds:
            if verbose: print "mark_txpts(): No CDS for %s." % gene_name
            continue
        # txpt_id = cds[gene_name].keys()[0]
        map_using_utr(flocs_map, gene_name, utrs, exon, cds, other,
                      sequences, rc_sequences)
    for gene_name in flocs_map:
        flocs_map[gene_name].check_coherency()
    return flocs_map


def flip_minus_strand_features(
            sequences, rc_sequences, chr_len, peaks,
            utr, exon, cds, other, no_peaks=False):
    for _feat in [utr, exon, cds, other]:
        for gene_name in _feat:
            for txpt_id in _feat[gene_name]:
                try:
                    if _feat[gene_name][txpt_id][0][6] == '+':
                        continue
                    chrm = _feat[gene_name][txpt_id][0][0]
                except: continue
                for _row in _feat[gene_name][txpt_id]:
                    old_right = int(_row[4])
                    _row[4] = chr_len[chrm] - (int(_row[3]) - 0) + 1
                    _row[3] = chr_len[chrm] - (old_right - 0)
    if no_peaks:
        return None
    for index, row in peaks.iterrows():
        if row['strand'] == '+': continue
        _new_right = chr_len[row['chrm']] - row['left'] + 1
        peaks.loc[index, 'left'] = chr_len[row['chrm']] - peaks.loc[index, 'right']
        peaks.loc[index, 'right'] = _new_right
    return peaks


def lack_parts(gene_name, txpt_id, exon):
    if gene_name not in exon:
        if verbose: print "mark_txpts(): Gene name %s not in exon dict." % gene_name
        return True
    if txpt_id not in exon[gene_name]:
        if verbose: print "mark_txpts(): Txpt id %s not in exon dict for gene %s." % (
            txpt_id, gene_name)
        return True
    if len(exon[gene_name][txpt_id]) == 0:
        if verbose: print "mark_txpts(): No exons."
        return True
    return False


def set_utr_spans(sorted_cds, utrs, cds_span, gene_name, txpt_id):
    # Set the UTR span as the end coordinates of the CDS,
    # or longer, if there is a UTR that extends further.
    left_utr_span = [int(sorted_cds[0][3]), int(sorted_cds[0][3])-1]
    right_utr_span = [int(sorted_cds[-1][4]), int(sorted_cds[-1][4])]
    given_left_utr = False
    given_right_utr = False
    if (gene_name in utrs) and (txpt_id in utrs[gene_name]) and (len(utrs[gene_name][txpt_id]) > 0):
        left_utr_iv = sorted(utrs[gene_name][txpt_id], key=lambda x: int(x[3]))
        right_utr_iv = sorted(utrs[gene_name][txpt_id], key=lambda x: int(x[4]))
        if int(left_utr_iv[0][3]) < int(cds_span[0]):
            left_utr_span = [int(left_utr_iv[0][3]), int(cds_span[0])]
            given_left_utr = True
        if int(right_utr_iv[-1][4]) > int(cds_span[1]):
            right_utr_span = [int(cds_span[1]), int(right_utr_iv[-1][4])]
            given_right_utr = True
    return (left_utr_span, right_utr_span, given_left_utr, given_right_utr)
    

def map_using_utr(flocs_map, gene_name, utrs, exon, cds, other,
                  sequences, rc_sequences):
    txpt_id = exon[gene_name].keys()[0]
    if lack_parts(gene_name, txpt_id, exon):
        return
    # Determine the UTR and CDS bounds.
    sorted_cds = sorted(cds[gene_name][txpt_id], key=lambda x: int(x[3]))
    cds_span = [int(sorted_cds[0][3]), int(sorted_cds[-1][4])]
    (left_utr_span, right_utr_span, given_left_utr, given_right_utr
     ) = set_utr_spans(sorted_cds, utrs, cds_span, gene_name, txpt_id)
#    sorted_exons = sorted(exon[gene_name][txpt_id], key=lambda x: int(x[3]))
    # Create flocs object for this gene and have it set
    # coordinate values relative to the pre-mRNA start, so the - strand can be used
    strand = sorted_cds[0][6]
    chrom = sorted_cds[0][0]
    exon_dict = {}
    exon_borders_in_seq = {}
    exon_num = 1
    # Start with the left utr region, if it exists.
    adjust = 0 if strand == '+' else 0
    left_utr_iv = [sorted_cds[0][0], left_utr_span[0]-1, left_utr_span[1]+adjust-1, '+']
    seq = seq_from_iv(left_utr_iv, sequences if strand == '+' else rc_sequences)
    exon_dict[exon_num] = left_utr_iv
    exon_borders_in_seq[exon_num] = [0, len(seq)]
    for ex in sorted_cds:
        slice_coord = [ex[0], ex[3] - 1, ex[4], '+']
        init_seq_len = len(seq)
        exon_dict[exon_num] = slice_coord
        #exon_list.append(slice_coord)
        seq += seq_from_iv(slice_coord, sequences if strand == '+' else rc_sequences)
        exon_borders_in_seq[exon_num] = [init_seq_len, len(seq)]
        if verbose: print 'mark_txpts(): set exon number %i as %s.' % (
            exon_num, str(slice_coord))
        exon_num += 1
    # Add right utr.
    init_seq_len = len(seq)
    right_utr_iv = [sorted_cds[0][0], right_utr_span[0]-1, right_utr_span[1], '+']
    seq += seq_from_iv(right_utr_iv, sequences if strand == '+' else rc_sequences)
    exon_dict[exon_num] = right_utr_iv
    exon_borders_in_seq[exon_num] = [init_seq_len, len(seq)]
    #except: print "Some oddity with exons."
    txpt_span = [left_utr_span[0]-1, right_utr_span[1]]
    if given_left_utr:
        txpt_span[0] = left_utr_span[0]
    if given_right_utr:
        txpt_span[1] = right_utr_span[1]
    # Get start, end, left_txpt_start and right_txpt_end.
    # Also load some sequences.
    # .gtf files are 1-based. Sequence slices are 0-based.
    # The left border of slices are (of course) included.
    # The right border of some exons x, say the third A in TAA, should be found at sequences[x-1:x].
    # To say the same thing, a stop codon in the .gtf at 647-649 will be found at sequences[646:649].
    # To find a start codon on the - strand located at 2167-2169 in a .gtf, the slice is sequences[2166:2169],
    # or, in other words, nt x:x+n (not including the final nt) on either strand is [x-1:x-1+n].
    # The downstream iv should start with the nt just after the stop, and be (say) 200 nt long.
    # In the .gtf, this would be stop_right_border+1:stop_right_border+1+200.
    # As a slice, this would therefore be stop_right_border:stop_right_border+200.
    # Upstream is start_left_border-200:start_left_border in the .gtf,
    # so start_left_border-200-1:start_left_border-1 as a slice.
#            upstream_iv = (exons[0][0], left_utr[0], left_utr[1] strand)  # Converted 1- to 0-based.
#            downstream_iv = (exons[0][0], right_utr[0], right_utr[1], strand)  # Converted 1- to 0-based.
    seq = extend_to_stop(cds_span, txpt_span, seq, sorted_cds, strand, sequences,
                         rc_sequences,
               exon_borders_in_seq, gene_name)
    flocs_map[gene_name] = flocs(
        txpt_id, chrom, txpt_span[0], txpt_span[1],
        cds_span[0], cds_span[1], strand, seq, exon_dict, gene_name,
        exon_borders_in_seq, given_left_utr, given_right_utr)
    flocs_map[gene_name].given_left_utr = given_left_utr
    flocs_map[gene_name].given_right_utr = given_right_utr
    #if verbose: print "mark_txpts(): Created flocs object %s" % str(flocs_map[gene_name])


def extend_to_stop(cds_span, txpt_span, seq, sorted_cds, strand, sequences,
                   rc_sequences,
                   exon_borders_in_seq, gene_name):
    if (cds_span[1] == txpt_span[1]) and (seq[-3:] not in ['TAA', 'TAG', 'TGA']):
        if verbose: print "Extending %s to get stop..." % gene_name
        _add = seq_from_iv(
            (sorted_cds[0][0], cds_span[1], cds_span[1] + 3, '+'),
            sequences if strand == '+' else rc_sequences)
        if _add in ['TAA', 'TAG', 'TGA']:
            txpt_span[1] += 3
            seq = seq + _add
            if verbose: print "Did it work? Expected stop: %s" % seq[-3:]
            try:
                exon_borders_in_seq[len(exon_borders_in_seq)][1] = len(seq)
            except:
                loggern.warn(
                    "%s: Missing an exon? %s" % (gene_name, str(exon_borders_in_seq)))
        else:
            if verbose: print "Did it work? Expected stop: %s" % _add
    return seq


def seq_from_iv(iv, sequences):
    """Returns a slice from 0-based, slice coordinates.
    """
    if iv[2] > len(sequences[iv[0]]):
        end = len(sequences[iv[0]])
    else:
        end = iv[2]
    a_seq = sequences[iv[0]][iv[1]:end]
    if iv[3] == "-":
        a_seq = rc(a_seq)
    return a_seq


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s


def get_sequences():
    fasta_filename = '/scratch/indexes/WS235.fa'
    sequences = dict(
        (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
    rc_sequences = dict(
        (p.name.split(' ')[0], rc(p.seq)) for p in HTSeq.FastaReader(fasta_filename))
    chr_lens = dict(
        [(name, len(sequences[name])) for name in sequences])
    return (sequences, rc_sequences, chr_lens)


def has_intron(txpt_id, gtf_df):
    if len(gtf_df[(gtf_df['transcript_id']==txpt_id) & (gtf_df['2']=='CDS')]) > 1:
        return True
    return False


def add_peak_locations_to_transcripts(peaks_l, txpts, chr_lens):
    all_added_peak_locs = []
    for gene_name in peaks_l:
        if gene_name in txpts:
            txpts[gene_name].find_peaks(peaks_l, chr_lens)
            all_added_peak_locs.append(len(txpts[gene_name].peak_locs))
        else:
            # print "Did not find {g} in txpts.".format(g=gene_name, )
            all_added_peak_locs.append(0)
    logger.info('''For genes with at least one peak:
    they have mean {mn} and median {md} peaks.'''.format(
        mn=np.mean([x for x in all_added_peak_locs if x>0]),
        md=np.median([x for x in all_added_peak_locs if x>0])
    ))
    if np.mean(all_added_peak_locs) == 0:
        logger.error("No peaks were assigned to targetted genes!")
    return txpts


def make_cdna_gtf(in_name, out_name):
    outf = open(out_name, 'w')
    with open(in_name, 'r') as f:
        header = next(f)
        outf.write(header)
        for li in f:
            if re.search('protein_coding', li) is not None:
                outf.write(li)


def load_signal_from_bedgraph(
        peaks, txpts, chr_lens,
        load_existing_data=False, do_combine_bedgraphs=False,
        bedgraphs_folder='data/wigs_coverage/',
        lib=None):
    if lib is not None:
        bedgraphs_folder = lib['bedgraphs_folder']
        txpts_obj_path = lib['txpt_obj_path']
    else:
        txpts_obj_path = 'data/for_heatmap.p'
    if load_existing_data:
        logger.info('%s: Loading data/for_heatmap.p holding the txpts object.' % (
            datetime.datetime.now().strftime('%Hh%Mm')))
        with open(txpts_obj_path, 'r') as f:
            txpts = pickle.load(f)
    else:
        logger.info('%s: Building a genomic array from %s.' % (
            datetime.datetime.now().strftime('%Hh%Mm'), bedgraphs_folder))
        ga = get_bedgraph(do_combine_bedgraphs=do_combine_bedgraphs,
                          bedgraphs_folder=bedgraphs_folder, lib=lib)
        logger.info('%s: Adding raw reads to txpts.' % (
            datetime.datetime.now().strftime('%Hh%Mm')))
        for gene_name in set(peaks['gene_name']):
            if gene_name not in txpts: continue
            txpts[gene_name].add_raw_reads_to_utr(ga, chr_lens)
            txpts[gene_name].add_raw_reads_to_all_peak_regions(ga)
        logger.info('%s: Saving txpts object to data/for_heatmap.' % (
            datetime.datetime.now().strftime('%Hh%Mm')))
        with open(txpts_obj_path, 'w') as f:
            pickle.dump(txpts, f)
    logger.info('Number of peaks: %i. Number of transcripts: %i.' % (
        len(peaks.index), len(txpts)
    ))
    locs = []
    for g in txpts:
        locs.append(len(txpts[g].peak_locs))
    logger.info("after loading, average number of peak locs {mn}".format(mn=np.mean(locs)))
    return peaks, txpts


def get_bedgraph(
        do_combine_bedgraphs=False, bedgraphs_folder='data/wigs/',
        lib=None):
    if lib is not None:
        bedgraphs_folder = lib['coverage_wigs']
        bedgraph_exp_plus = lib['bedgraph_exp_plus']
        bedgraph_exp_minus = lib['bedgraph_exp_minus']
    else:
        bedgraph_exp_plus = bedgraphs_folder + 'both_fbfs_plus.bed'
        bedgraph_exp_minus = bedgraphs_folder + 'both_fbfs_minus.bed'
    bedgraphs_folder = bedgraphs_folder.rstrip('/') + '/'
    if do_combine_bedgraphs: combine_bedgraphs(bedgraphs_folder=bedgraphs_folder)
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    with open(bedgraph_exp_plus, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] = float(s[3])
    with open(bedgraph_exp_minus, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] = float(s[3])
    return ga


def combine_bedgraphs(bedgraphs_folder='data/wigs_five_prime/'):
    ga = {}
    bedgraphs_folder = bedgraphs_folder.rstrip('/') + '/'
    for filename_list in [
            (bedgraphs_folder + 'fbf1_reads_plus.bed',
             bedgraphs_folder + 'fbf1_reads_minus.bed',
             'combined_fbf1.txt'),
            (bedgraphs_folder + 'fbf2_reads_plus.bed',
             bedgraphs_folder + 'fbf2_reads_minus.bed',
             'combined_fbf2.txt')]:
        peaks_filename = filename_list[2]
        scatterplot_correlation_by_wig.load_bedgraph(filename_list, ga, use_key=peaks_filename)
    ga['combined'] = HTSeq.GenomicArray(chroms='auto', stranded=True)
    for iv, score in ga['combined_fbf1.txt'].steps():
        ga['combined'][iv] += score
    for iv, score in ga['combined_fbf2.txt'].steps():
        ga['combined'][iv] += score
#    with open('temp_ga.p', 'w') as f:
#        pickle.dump(ga, f)
    ga['combined'].write_bedgraph_file(bedgraphs_folder + 'both_fbfs_plus.bed', '+')
    ga['combined'].write_bedgraph_file(bedgraphs_folder + 'both_fbfs_minus.bed', '-')
    return ga


def get_args():
    parser = argparse.ArgumentParser(description='''
        Determine the locations of features in genes across the genome.''')
    parser.add_argument('-i', '--input',
                        default='new_pypeaks_fdr1_negip_local_five_reps/combined_fbf.txt',
                        help='''Input peaks file (Required).''')
    parser.add_argument('-l', '--load_txpts',
                        help='Load marked transcript information.',
                        default=False, action='store_true')
    parser.add_argument('-m', '--heat_map',
                        default=False, action='store_true',
                        help='Just make the heatmap.')
    parser.add_argument('-g', '--load_gtf',
                        help='Load objects created by parsing a gtf file. Implies load_txpts is False.',
                        default=False, action='store_true')
    parser.add_argument('-c', '--config',
                        help='Directory holding config.py',
                        default='None')
#    parser.add_argument('-w', '--continuous_signal',
#                        default=False, action='store_true',
#                        help='Load true mapping data to display peak signal (slow).')
    args = parser.parse_args()
    return args


def get_gtf(args, lib=None):
    if lib is None:
        gtfname = '../clip/lib/gtf_one_txpt_per_gene.txt'
        gtf_pickle = 'data/feat_loc/gtf_data.p'
    else:
        gtfname = lib['gtf_one_txpt_per_gene']
        gtf_pickle = lib['gtf_pickle']
    if not args.load_gtf:  # Not loading from a pickle object.
        #if not os.path.exists(gtfname):        
        (utr, exon, cds, other) = load_gtf(
                gtfname, return_split_lines=True)
        with open(gtf_pickle, 'w') as f:
            pickle.dump((utr, exon, cds, other), f)
    else:
        with open(gtf_pickle, 'r') as f:
            (utr, exon, cds, other) = pickle.load(f)
    return (utr, exon, cds, other)

def create_one_txpt_per_gene(gtf_in, gtf_out):
    gtf = pandas.read_csv(gtf_in, sep='\t')
    gene_to_txpt = zip(gtf.gene_name, gtf.transcript_id)
    #txpt_to_left = zip(gtf.transcript_id, gtf.3)
    #txpt_to_right =zip(gtf.transcript_id, gtf.4)
    min_left = collections.defaultdict(list)
#    for t in gtf.transcript_id:
        

def get_peaks(filename, sequences, rc_sequences, chr_lens,
                      utr, exon, cds, other):
    peaks = pandas.read_csv(filename, sep='\t')
    logger.info("get_peaks(): Loaded %i peaks in %i genes from %s." % (
        len(peaks['gene_name'].tolist()),
        len(list(set(peaks['gene_name'].tolist()))), filename))
    peaks = flip_minus_strand_features(
        sequences, rc_sequences, chr_lens, peaks, utr, exon, cds, other)
    # - strand features are still listed on the - strand,
    # but their genomic coordinates are relative to the 3' end
    # of the + strand now. They can be treated as + strand,
    # except that sequences should be read from rc_sequences.
    return peaks


def get_txpts(args, sequences, rc_sequences, peaks, utr, exon, cds, other,
              lib=None):
    print 'get_txpts'
    if lib is None:
        feat_loc = 'data/feat_loc'
    else:
        feat_loc = lib['feat_loc_dir']
    if not os.path.exists(feat_loc):
        os.system('mkdir ' + feat_loc)
    if not args.load_txpts:
        print 'marking txpts...'
        txpts = mark_txpts(sequences, rc_sequences, peaks, utr, exon, cds, other)
        print len(txpts)
        print 'get_txpts after mark_txpts'
        print lib
        for t in txpts:
            txpts[t].check_coherency()
            txpts[t].motif = lib['motif']
            # flocs_map[gene_name].find_features()
            print "finding features"
            txpts[t].find_features()
        #add_peak_locations_to_transcripts(peaks, txpts, chr_lens)
        with open(feat_loc + '/%s.p' % os.path.basename(args.input).partition('.txt')[0], 'w') as f:
            pickle.dump(txpts, f)
    else:
        with open(feat_loc + '/%s.p' % os.path.basename(args.input).partition('.txt')[0], 'r') as f:
            txpts = pickle.load(f)
    logger.info("\tLoaded %i transcripts with features." % (len(txpts)))
    return txpts


def get_data(args, lib=None):
    logger.info('%s: Reading sequence information.' % datetime.datetime.now().strftime('%Hh%Mm'))
    print "Getting sequences..."
    (sequences, rc_sequences, chr_lens) = get_sequences()
    logger.info('%s: Getting GTF.' % datetime.datetime.now().strftime('%Hh%Mm'))
    print "Getting GTF..."
    (utr, exon, cds, other) = get_gtf(args, lib=lib)
    logger.info('%s: Getting peaks from %s.' % (
        datetime.datetime.now().strftime('%Hh%Mm'), args.input))
    peaks = get_peaks(args.input, sequences, rc_sequences, chr_lens,
                      utr, exon, cds, other)
    print '---- getting txpts...'
    logger.info('%s: Getting transcripts.' % datetime.datetime.now().strftime('%Hh%Mm'))
    txpts = get_txpts(
        args, sequences, rc_sequences, peaks, utr, exon, cds, other, lib=lib
    )
    print '---- finished getting txpts.'
    peaks_d = to_dict(peaks)
    logger.info('%s: Adding peak locations.' % datetime.datetime.now().strftime('%Hh%Mm'))
    txpts = add_peak_locations_to_transcripts(peaks_d, txpts, chr_lens)
    locs = []
    for g in txpts:
        locs.append(len(txpts[g].peak_locs))
    logger.info("locs {mn}".format(mn=np.mean(locs)))
    peaks_l = []
    for gene_name in peaks_d:
        peaks_l.extend(peaks_d[gene_name])
    peaks = pandas.DataFrame(peaks_l)
    logger.info('%s: Loading signal from bedgraph.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    peaks, txpts = load_signal_from_bedgraph(
        peaks, txpts, chr_lens, load_existing_data=args.load_txpts,
        lib=lib)
    return peaks, txpts, chr_lens


def to_dict(df):
    _li = {}
    for index, row in df.iterrows():
        _li.setdefault(row['gene_name'], [])
        _li[row['gene_name']].append(row.to_dict())
    return _li


def make_figs(peaks, txpts, chr_lens, args, input_dir, lib=None):
    locs = []
    #pandas.DataFrame.sort(columns='height', ascending=False, inplace=True)
    peaks.sort(columns='height', ascending=False, inplace=True)
    output_dirname = 'figs/%s/' % os.path.dirname(args.input)
    for g in txpts:
        locs.append(len(txpts[g].peak_locs))
    logger.info("Assigned %i peaks and %i transcripts." % (len(peaks.index), len(txpts)))
#    logger.info("before loading, average number of peak locs {mn}".format(mn=np.mean(locs)))
#    args.load_txpts = True
    print "Assigned %i peaks and %i transcripts." % (len(peaks.index), len(txpts))
    logger.info('%s: Peak in gene line plot.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    peak_in_gene_line_plot.normalize_distances(txpts)
    logger.info('%s: Getting average positions.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    (ave_peaks, ave_fbes, ave_negs, ave_polyA, ave_peak_v_polyA, ave_fbe_v_polyA,
     ave_highest_peak, ave_secondary_peaks) = peak_in_gene_line_plot.get_ave_pos(txpts)
    #logger.info('%s: Peak vs FBE line plot.' % (
    #    datetime.datetime.now().strftime('%Hh%Mm')))
    #peak_in_gene_line_plot.peak_vs_fbe(txpts)
    #logger.info('%s: Features in normalized gene line plot.' % (
    #    datetime.datetime.now().strftime('%Hh%Mm')))
    heatmap_of_raw_signal(peaks, txpts, output_dirname=output_dirname,
                          include_motif=True)
    peak_in_gene_line_plot.plot_features(
       (ave_peaks, ave_fbes, ave_negs, ave_polyA, ave_highest_peak, ave_secondary_peaks),
       output_filename='figs/%s/features_in_normalized_gene.pdf' % os.path.dirname(input_dir))
    #logger.info('%s: Plotting UTR line plot.' % (
    #    datetime.datetime.now().strftime('%Hh%Mm')))
    peak_in_gene_line_plot.plot_utr(ave_peak_v_polyA, ave_fbe_v_polyA,
              output_filename='figs/%s/features_in_utr.pdf' % os.path.dirname(input_dir))
    logger.info('%s: Heatmap by gene length, only peak ranges.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    heatmap_by_gene_length_just_peak_ranges(
        peaks, txpts, output_dirname='figs/%s/' % os.path.dirname(input_dir))
    logger.info('\n***\n%s: Creating heatmaps.\n***' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    # logger.info('%s: Peak regions heatmaps.' % (
    #     datetime.datetime.now().strftime('%Hh%Mm')))
    heatmap_of_peak_region(peaks, txpts, output_dirname=output_dirname)
    # logger.info('%s: Intrapeak distances.' % (
    #     datetime.datetime.now().strftime('%Hh%Mm')))
    output_heatmap_figures.plot_intrapeak_distances(peaks, txpts, output_dirname=output_dirname)
    # logger.info('%s: Raw signal heatmaps.' % (
    #     datetime.datetime.now().strftime('%Hh%Mm')))

if __name__ == '__main__':
    args = get_args()
    logger = logging.getLogger(__name__)
    if args.config != 'None':
        sys.path.insert(0, args.config)
        import config
        lib = config.config()
    else:
        lib = {'gtf': 'lib/gtf_with_names_column.txt',
               'gtf_one_txpt_per_gene': 'lib/gtf_one_txtp_per_gene.txt'}
    # if args.load_gtf:
    #     args.load_txpts = True
    (sequences, rc_sequences, chr_lens) = get_sequences()
#    make_cdna_gtf('../clip/lib/gtf_with_biotype_column.txt',
#                  '../clip/lib/gtf_only_mrna.txt')
    input_dir = os.path.dirname(args.input) + '/' #'../clip/pypeaks_fdr1_negip_local/'
    if not args.load_gtf:
        create_file_of_longest_transcript_per_gene(
            in_filename=lib['gtf'],
            out_filename=lib['gtf_one_txpt_per_gene'])#'../clip/lib/gtf_one_txpt_per_gene.txt')
    sys.exit()
    if not os.path.exists('./data'): os.system('mkdir data')
    if not os.path.exists('./data/feat_loc'): os.system('mkdir data/feat_loc')
    if args.heat_map:
        output_filename = 'figs/%s/feat_heatmap.pdf' % os.path.dirname(args.input)
        print "Making heatmap in %s (after loading pickled txpt data)..." % output_filename
        with open('data/feat_loc/%s.p' % os.path.basename(args.input).partition('.txt')[0], 'r') as f:
            txpts = pickle.load(f)
        peaks = pandas.read_csv(args.input, sep='\t')
        if not os.path.exists('figs/%s/' % os.path.dirname(args.input)):
            os.system('mkdir ' + 'figs/%s/' % os.path.dirname(args.input))
        load_signal_from_bedgraph(
            peaks, txpts, chr_lens, load_existing_data=args.load_data)
        heatmap_by_gene_length_just_peak_ranges(
            peaks, txpts, output_dirname='figs/%s/' % os.path.dirname(args.input))
        sys.exit()
    (peaks, txpts, chr_lens) = get_data(args, lib=lib)
    make_figs(peaks, txpts, chr_lens, args, input_dir)
