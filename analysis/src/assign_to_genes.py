import HTSeq
import sys
import os
import re
import argparse
import collections
import glob
import time
import collections

import config

# Following:
# http://www.bioconductor.org/help/workflows/rnaseqGene/
#samfile = sys.argv[1]

# HTSeq segment from:
# http://www-huber.embl.de/users/anders/HTSeq/doc/counting.html#counting

skip ='''
def a(lib, infile):
    df = pandas.read_csv(gtf['gtf'], sep='\t')
    id_to_name = dict(zip(df.gene_id, df.gene_name))
    peaks = pandas.read_csv(infile, sep='\t')
    peaks['oo_rna_seq_gene'] = [
        name_to_counts['oo'][x] for x in peaks.gene_name]
    peaks['modencode_rna_seq_gene'] = [
        name_to_counts['modencode'][x] for x in peaks.gene_name]
    peaks['oo_rna_seq_in_peak'] = []
        '''

def get_gtf(lib, gtf_file=None):
    if gtf_file is not None:
        gtf_path = gtf_file
    else:
        gtf_path = lib['gtf_raw']
    print gtf_path
    gtf_file = HTSeq.GFF_Reader(gtf_path)
    features = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    for feature in gtf_file:
        if feature.type == 'exon':
            features[feature.iv] += feature.attr['gene_name']
    return features


def map_folder(bed_dir, features):
    for bed_filename in glob.glob(bed_dir + '/*.bed'):
        print "Assigning reads from %s..." % bed_filename
        map_bed(bed_filename, features)


def read_bed(bed_file):
    switch_strand = {'-': '+', '+': '-'}
    with open(bed_file, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            # num_reads = int(re.match('.*n=(\d+).*', li).groups()[0])
            yield (
                HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), s[5]), #switch_strand[s[5]]),
#                 HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), switch_strand[s[5]]),
                1)


def map_bed(bed_filename, features):
    start_time = time.clock()
#    counts = collections.collections.Counter()
    counts = collections.defaultdict(int)
    for line_num, (read_iv, number_of_reads), in enumerate(read_bed(bed_filename)):
        if not line_num % 100000:
            print line_num
            elapsed = time.clock() - start_time
            per_ten_million = 1e6 * elapsed/max([1, float(line_num)])
            print "Time elapsed: %f. Seconds per million reads %s." % (
                elapsed, per_ten_million)
#        if not almnt.aligned:
#            counts['_unmapped'] += 1
#            continue
        gene_ids = set()
        for iv, val in features[read_iv].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[gene_id] += number_of_reads
        elif len(gene_ids) == 0:
            counts['_no_feature'] += number_of_reads
        else:
            counts['_ambiguous'] += number_of_reads
    print 'sum %s' % bed_filename
    print sum([counts[x] for x in counts])
    try:
        print 'no feature %i ambiguous %i' % (counts['_no_feature'], counts['_ambiguous'])
    except:
        print 'either no ambiguous or no "no_feature" reads.'
    output_filename = 'counts/' + os.path.basename(bed_filename).partition('.bed')[0]
    output_filename += '_counts.txt'
    print output_filename
    if not os.path.exists('counts/'):
        os.system('mkdir counts')
#    if os.path.exists(output_filename):
#        continue
    with open(output_filename, 'w') as f:
        for gene_id in counts:
            f.write("{gid}\t{cts}\n".format(
                gid=gene_id, cts=counts[gene_id])
            )


def fill_in_gaps(folder_name, lib):
    # if gtf_file is None:
    #     gtf_path = gtf_file
    # else:
    #     gtf_path = lib['gtf_raw']
    wb_ids = gene_names_in_gtf(lib['gtf_raw'])
    counts = {}
    known_genes = set()
    for fname in glob.glob(folder_name + '/*'):
        print fname
        counts[fname] = collections.defaultdict(int)
        with open(fname, 'r') as f:
            for li in f:
                s = li.rstrip('\n').split('\t')
                known_genes.add(s[0])
                counts[fname][s[0]] += int(s[1])
    # for fname in counts:
    #     for gene in list(wb_ids - set(counts[fname].keys())):
    #         counts[fname][gene] = 0
    gene_output_order = list(known_genes)#counts[counts.keys()[0]].keys()
    for fname in counts:
        for gene in known_genes:
            if gene not in counts[fname]:
                counts[fname][gene] = 0
        with open(fname, 'w') as f:
            for gene in gene_output_order:
                f.write('{i}\t{n}\n'.format(i=gene, n=counts[fname][gene]))


def create_combined_file_of_counts(
        counts_folder, output_filename='combined_counts.txt'):
    import glob
    import collections
    all = collections.defaultdict(dict)
    s = {}
    genes = set()
    output_order = []
    for fname in glob.glob(counts_folder + '/*.txt'):
        f = open(fname, 'r')
        for line in f:
            s = line.rstrip('\n').split('\t')
            all[s[0]][os.path.basename(fname)] = s[1]
        output_order.append(os.path.basename(fname))
    li = 'gene\t' + '\t'.join(output_order) + '\n'
    for gene in all:
        li += '%s\t' % gene
        for k in output_order:
            li += '%s\t' % all[gene][k]
        li += '\n'
    with open(output_filename, 'w') as f:
        f.write(li)
        # s[fname] = dict([line.rstrip('\n').split('\t') for line in f])
        # genes.add(set(s.keys()))


def gene_names_in_gtf(gtf_fname):
    wb_ids = set()
    with open(gtf_fname, 'r') as f:
        for li in f:
            m = re.match('.*gene_name "([^"]+)";', li)
            if m is not None:
                wb_ids.add(m.groups()[0])
    return wb_ids


def wb_to_public_name(gtf_fname):
    public_names = {}
    with open(gtf_fname, 'r') as f:
        for li in f:
            m = re.match('.*gene_name "([^"]+)";', li)
            if m is not None:
                name_m = re.match('.*gene_name "([^"]+)";', li)
                if name_m is not None:
                    public_names[m.groups()[0]] = name_m.groups()[0]
    return public_names


def set_of_all_gene_names_in_folder(folder_name):
    gene_names = set()
    for fname in glob.glob(folder_name + '/*'):
        with open(fname, 'r') as f:
            for li in f:
                gene_names.add(li.rstrip('\n').split('\t')[0])
    return gene_names


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed', help="""Folder of bed files.""")
    parser.add_argument('-c', '--config', help="""config.ini""")
    args = parser.parse_args()
    return args


def run(args):
    lib = config.config(args.config)
    features = get_gtf(lib, gtf_file=lib['gtf_raw'])  #'./genomes/Saccharomyces_cerevisiae.EF4.70.gtf')
    map_folder(args.bed, features)
    fill_in_gaps('counts', lib)
    create_combined_file_of_counts(
        'counts', output_filename='combined_counts.txt')


if __name__ == '__main__':
    """Create a combined_counts.txt file of the reads/gene for a 
    folder of .bed files. Requires only one value from lib, lib['gtf_raw'].
    """
    run(get_args())
