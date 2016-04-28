import os
import glob
import re
import sys
import pandas
import HTSeq

def get_sequences(combined):
    fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for index, peak_row in combined.iterrows():
        start = combined.loc[index, 'left']
        end = combined.loc[index, 'right']
        chrm = combined.loc[index, 'chrm']
        seq = sequences[chrm][start:end]
        #print "%s:%i-%i: seq %s" % (chrm, start, end, seq)
        if combined.loc[index, 'strand'] == '-':
            seq = rc(seq)
        combined.loc[index, 'seq'] = seq
        
def complement(s):
    #basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)

def rc(s):
    s = s[::-1]
    s = complement(s)
    return s

def run_dreme(combined, label, top_level_dir):
    out_dir = 'data/' + os.path.basename(top_level_dir.rstrip('/'))
    out_dir += '/dCLIP_fasta'
    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)
        os.system('mkdir ' + out_dir + '/dreme_results')
    fasta_filename = '%s/%s.fa' % (out_dir, label)
    as_fasta = ""
    for index, peak_row in combined.iterrows():
        start = combined.loc[index, 'left']
        end = combined.loc[index, 'right']
        chrm = combined.loc[index, 'chrm']
        name = str(index)
        name += ":%s:%i:%i" % (chrm, start, end)
        as_fasta += ">%s\n%s\n" % (name, combined.loc[index, 'seq'])
    with open(fasta_filename, 'w') as f:
        f.write(as_fasta)
    dreme_out_dir = out_dir + '/dreme_results/%s' % label
    os.system('/home/dp/meme/bin/dreme -e 0.01 -p %s -norc -maxk 10 -oc %s' % (
        fasta_filename, dreme_out_dir))

def split_by_state_and_run_dreme(peaks, raw, top_level_dir):
    for index, row in peaks.iterrows():
        peak_id = peaks.loc[index, 'id']
        print 'peak id: %s' % str(peak_id)
        state = raw[raw['id']==peak_id]['state'].values[0]
        peaks.loc[index, 'state'] = state
    no_change = peaks[peaks['state']==1]
    fbf1_up = peaks[peaks['state']==2]
    fbf2_up = peaks[peaks['state']==0]
    run_dreme(fbf1_up, 'fbf1_enriched', top_level_dir)
    run_dreme(fbf2_up, 'fbf2_enriched', top_level_dir)
    run_dreme(no_change, 'no_difference', top_level_dir)


def split_and_assign_to_gene(peaks, gtf_df):
    peaks['gene_name'] = ''
    chrms = dict(peaks['chrm'].value_counts()).keys()
    for chrm in chrms:
        for strand in ['+', '-']:
            print "chrm %s, strand %s" % (chrm, strand)
            assign_to_gene(peaks, chrm, strand, gtf_df)


# Find the closet gene for each peak, and resolve ties.   
def assign_to_gene(peaks, chrm, strand, gtf_df, verbose=False):
    peaks_sub = peaks[(peaks['chrm']==chrm) & (peaks['strand']==strand)]
    sub = gtf_df[(gtf_df['0']==chrm) & (gtf_df['6']==strand)]
    sub = sub[sub['2']=='transcript']
    for index, apeak in peaks_sub.iterrows():
        chrm = peaks.loc[index, 'chrm']
        strand = peaks.loc[index, 'strand']
        if not index % 100:
            print "Assigning gene %i/%i for chrm/strand %s/%s" % (
                index, len(peaks['strand']), chrm, strand)
        #index_left = bisect.bisect_left(sub_list_left, apeak.left - 1e3)
        #index_right = bisect.bisect_left(sub_list_right, apeak.left + 1e3)
        asub = sub[(abs(sub['3'] - apeak.left) < 1e5
                    ) | (abs(sub['4'] - apeak.left) < 1e5)]
        genes_for_peak = find_nearby_genes(asub, apeak)
        if verbose:
            print "assign_to_gene(): .genes_for_peak="
            for a_gene in genes_for_peak:
                print "\n%s\n" % str(a_gene)
        closest = (1e4, None)
        ties = []
        for _gene in genes_for_peak:
            dist = get_distance_pandas(_gene, apeak)
            if verbose:
                print "\n***_gene %s: dist: %i***\n" % (_gene.gene_name, dist) 
            if abs(dist) < abs(closest[0]):
                ties = []
                closest = (dist, _gene)
                if verbose:
                    print "Set closest to %s, dist %i" % (_gene.gene_name, dist)
            elif abs(dist) == abs(closest[0]):
                ties.append((dist, dict(_gene)))
                ties.append((dist, dict(closest[1])))
        if len(ties) > 0:
            for tup in ties:
                dist = tup[0]
                _gene = tup[1]
                #dist = get_distance_pandas(_gene, apeak)
                if dist == closest[0]:
                    closest = resolve_tie(ties)
                    break
        peaks.loc[index, 'gene_dist'] = closest[0]
        if closest[1] is not None:
            peaks.loc[index, 'gene_name'] = closest[1]['gene_name']
        else:
            peaks.loc[index, 'gene_name'] = 'Unknown'


def find_nearby_genes(sub, apeak):
    """Assigns a gene for a given peak.
    """
    cutoff = 1e3
    genes_for_peak = []
    for index, _gene in sub.iterrows():
        dist = get_distance_pandas(_gene, apeak)
        if abs(dist) < cutoff:
            genes_for_peak.append(_gene)
    return genes_for_peak

def get_distance_pandas(_gene, apeak):
    """Used to assign peaks to a gene.
    """
    if _gene['4'] < apeak.left:
        return apeak.left - _gene['4']
    if apeak.right < _gene['3']:
        return apeak.right - _gene['3']
    return 0

def resolve_tie(ties):
    biotypes = []
    for tup in ties:
        _gene = tup[1]
        m = re.search('gene_biotype "([^"])"', _gene['8'])
        if m is not None:
            biotypes.append((tup[0], tup[1], m.group(1)))
        else:
            biotypes.append((tup[0], tup[1], 'Unknown'))
    # If there is a non-coding RNA, assign to that.
    non_coding = []
    for tup in biotypes:
        if tup[2] != 'protein_coding':
            non_coding.append(tup)
    if len(list(non_coding)) == 1:
        return (non_coding[0][0], non_coding[0][1])
    if len(list(non_coding)) > 1:
        # Multiple non-coding RNAs. Pick randomly.
        return (non_coding[0][0], non_coding[0][1])
    if len(list(non_coding)) == 0:
        # No non-coding RNA. Pick randomly.
        return (ties[0][0], ties[0][1])

def add_state(peaks, raw):
    for index, row in peaks.iterrows():
        peak_id = peaks.loc[index, 'id']
        state = raw[raw['id']==peak_id]['state'].values[0]
        peaks.loc[index, 'state'] = state

def write_differential_regions_with_genes(peaks):
    no_change = peaks[peaks['state']==1]
    fbf1_up = peaks[peaks['state']==2]
    fbf2_up = peaks[peaks['state']==0]
    no_change.to_csv('data/dCLIP_no_change.txt', sep='\t')
    fbf1_up.to_csv('data/dCLIP_fbf1_up.txt', sep='\t')
    fbf2_up.to_csv('data/dCLIP_fbf2_up.txt', sep='\t')

if __name__ == '__main__':
    top_level_dir = sys.argv[1]
    #top_level_dir = 'pypeaks_fdr1_negip_local/'
    in_dir = 'data/' + os.path.basename(top_level_dir.rstrip('/'))
    summary_file = in_dir + '/dCLIP_out/dCLIP_summary.bed'
    with open(summary_file, 'r') as f:
        next(f)
        out_li = "chrm\tleft\tright\tid\tstrength\tstrand\tleftx\trightx\tunk\n"
        for li in f:
            out_li += li
    with open(summary_file.rstrip('.bed') + '_new_header', 'w') as f:
        f.write(out_li)
    # The dCLIP_output.txt file contains the state information.
    raw_dCLIP_file = in_dir + '/dCLIP_out/dCLIP_output.txt'
    dCLIP_output = pandas.read_csv(raw_dCLIP_file, sep='\t')
    new_in = summary_file.rstrip('.bed') + '_new_header'
    peaks = pandas.read_csv(new_in, sep='\t')
    get_sequences(peaks)
    split_by_state_and_run_dreme(peaks, dCLIP_output, top_level_dir)
    gtf_df = pandas.read_csv('lib/gtf_with_names_column.txt', sep='\t')
    add_state(peaks, dCLIP_output)
    split_and_assign_to_gene(peaks, gtf_df)
    write_differential_regions_with_genes(peaks)
