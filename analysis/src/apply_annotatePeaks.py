import glob
import pandas
import sys
import HTSeq
import os
import argparse
import collections
import re
sys.path.insert(0, '/groups/Kimble/Aman Prasad/clip2/src/')
from annotatedPeaks import annotatedPeaks
import config

del sys.path[0]

def add_info(fname, gtf, sequences, ga, gtf_d):
    p = annotatedPeaks(file=fname)
    p.add_biotype(gtf)
    p.add_seqs(sequences)
    p.add_fbe()
    p.copy_columns()
    p.add_location_from_integral(gtf_d, ga)
    p.write(fname)


def load_bedgraph(fname):
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    plus_file = fname.partition('.wig')[0] + '_+.wig'
    add_strand_to_ga_from_bedgraph_file(plus_file, ga, '+')
    minus_file = fname.partition('.wig')[0] + '_-.wig'
    add_strand_to_ga_from_bedgraph_file(minus_file, ga, '-')
    return ga


def add_strand_to_ga_from_bedgraph_file(fname, ga, strand):
    with open(fname, 'r') as f:
        next(f)
        for line in f:
            #if line[0] != 'M': continue
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), strand)] = float(s[3])
    return ga


def get_val(ga, iv):
    return np.max(np.fromiter(ga[iv], dtype=np.float))


def get_bed_size(fname):
    return float(len(open(fname).readlines()))


def create_raw_reads_columns(lib, rootdir):
    sizes = get_bed_sizes(lib)
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            path = os.path.join(subdir, file)
            p = annotatedPeaks(file=path)
            p.absolute_read_number(sizes)
            p.write(path)

def make_ratio_column_from_raw_reads(lib, rootdir):
    sizes = get_bed_sizes(lib)
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            path = os.path.join(subdir, file)
            print path
            p = annotatedPeaks(file=path)
            p.absolute_read_number(sizes)
            p.set_sum(
                to_sum=['unnorm_reads_fbf1_rep_1', 'unnorm_reads_fbf1_rep_2', 'unnorm_reads_fbf1_rep_3'],
                        summed_col='unnorm_reads_fbf1_reads')
            p.set_sum(
                to_sum=['unnorm_reads_fbf2_rep_1', 'unnorm_reads_fbf2_rep_2', 'unnorm_reads_fbf2_rep_3'],
                        summed_col='unnorm_reads_fbf2_reads')
            if re.search('fbf1', path):
                p.set_ratio(col1='unnorm_reads_fbf1_reads', col2='unnorm_reads_fbf1_n2')
            if re.search('fbf2', path):
                p.set_ratio(col1='unnorm_reads_fbf2_reads', col2='unnorm_reads_fbf2_n2')
            p.write(path)
            print p


def get_bed_sizes(lib):
    reps_labels = ['fbf1_rep_1', 'fbf1_rep_2', 'fbf1_rep_3',
        'fbf2_rep_1', 'fbf2_rep_2', 'fbf2_rep_3',
            'rna_seq', 'fbf1_n2', 'fbf2_n2']
    size = collections.defaultdict(int)
    #  fbf1_rep_1: %(top)s/fbf1/bedgraph_norm/exp_fbf1_TGGC
    #  read_beds: %(top)s/bed_collapsed/
    for bed in [x for x in lib if re.search('_bed$', x)]:
        if type(lib[bed]) != type(''):
            continue
        label = bed.partition('_bed')[0]
        path = lib[bed] + '.bed'
        print "***" + path
        if os.path.exists(path):
            size[label] = get_bed_size(path)
        else:
            print "Expected file %s, but could not find..." % path
    print size
    return size


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-r', '--ratio_is_raw_reads',
                        default=False, action='store_true')
    parser.add_argument('-c', '--config_ini')
    args = parser.parse_args()
    lib = config.config(args.config_ini)
    rootdir = args.input
    gtf = pandas.read_csv('lib/gtf_with_names_column.txt', sep='\t')
    fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    gtf_r = gtf.to_dict('records')
    gtf_d = collections.defaultdict(list)
    for row in gtf_r:
        gtf_d[row['gene_name']].append(row)
    #gtf_one_per = dict([(x['gene_name'], x) for x in gtf_one_per])
    ga = load_bedgraph('bedgraph_norm/combined_fbf1.wig')
    for subdir, dirs, files in os.walk(rootdir):
        print "*" * 14
        print subdir
        print "---"
        print dirs
        print "---"
        print files
        for file in files:
            print os.path.join(subdir, file)
            add_info(os.path.join(subdir, file),
                     gtf, sequences, ga, gtf_d)
    if args.ratio_is_raw_reads:
        make_ratio_column_from_raw_reads(lib, rootdir)
