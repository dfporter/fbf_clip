"""
Load txpt objects and read locations within txpts.

For each txpt: call peaks by read shuffling (in negative and positive) and by CWT.

Rank by p value for the shuffle.

Apply a BH FDR.

INPUT:
Config file.

OUTPUT:
Outputs to lib['clusters_dir'] from the config.ini file.

Run as;
python find_peaks_by_permutations.py -c <directory with config.py> --use lib

From config.ini, this uses:
lib[control_bed1]...
lib[expt_bed1]...
lib['gtf_raw']
lib['read_beds']
lib['gtf_with_names']
exp_bedgraph_minus
exp_bedgraph_plus

cluster_tables_to_peak_files:
gtf_with_names
gtf_raw
clusters_dir
permutation_peaks_dir

subpeaks_ui.py
lib['coverage_wigs.]
exp_bedgraph_minus
exp_bedgraph_plus
lib['chr_sizes']

determine_feature_locations
gtf_pickle
feat_loc_dir
gtf_one_txpt_per_gene
exp_bedgraph_minus
exp_bedgraph_plus
coverage_wigs = bedgraphs_folder
txpt_obj_path
bedgraphs_folder
"""

import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys
import numpy as np
import pandas
import collections
import re

import determine_feature_locations
import subpeaks
import put_reads_in_gene
import cluster_reads
import p_values_of_clusters


def read_args():
    parser = argparse.ArgumentParser(description='''
        Determine the locations of features in genes across the genome.''')
    parser.add_argument(
        '-n', '--no_user_input', action='store_true', default=False,
        help='Do not request user input or run as a loop (Default: False).'
    )
    parser.add_argument('-i', '--input',
                        default='dummy_value',
                        help='''Don't set this value at present.''')
    parser.add_argument('-c', '--config_dir',
                        default='config/',
                        help='''Directory holding a config.py file.''')

    parser.add_argument('-v', '--no_ui',
                        help='No user input.',
                        default=False, action='store_true')
    parser.add_argument('-l', '--load_txpts',
                        help='Load marked transcript information.',
                        default=False, action='store_true')
#    parser.add_argument('-m', '--heat_map',
#                        default=False, action='store_true',
#                        help='Just make the heatmap.')
    parser.add_argument('-g', '--load_gtf',
                        help='Load objects created by parsing a gtf file. Implies load_txpts is False.',
                        default=False, action='store_true')
    parser.add_argument('-u', '--use',
                        help="""
Which fbf to use (default: both). Can also be lib (if not
analyzing FBF, in which case lib[exp_bedgraph_plus] and
lib[exp_bedgraph_mins] are used) or fbf1 or fbf2.""",
                        default='both')
    args = parser.parse_args()
    return args


def get_bed_files(bed_folder='beds/', gtf=None, args=None, lib=None):
    ht_exons = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    gtfd = collections.defaultdict(set)
    for feature in gtf:
        if feature.type == 'exon':
            ht_exons[feature.iv] += feature.name
            gtfd[feature.name].add(feature)
    control_names = [lib[x] for x in lib.keys() if re.match('control_bed.*', x) is not None]
    exp_names = [lib[x] for x in lib.keys() if re.match('exp_bed.*', x) is not None]
    print control_names
    print exp_names
    all_names = control_names + exp_names
    bed_file_list = [
        "/".join([bed_folder, x]) for x in control_names + exp_names
        # [
        #     lib['control_bed1'], lib['control_bed2'], lib['control_bed3'],
        # lib['control_bed4'], lib['exp_bed1'], lib['exp_bed2'],
        # lib['exp_bed3'], lib['exp_bed4']]
    ]
    reads_by_gene = {}
    counts_by_gene = {}
    ga = {}
    for fname in bed_file_list:
        (_ga, reads_by_gene[fname], counts_by_gene[fname]) = read_a_bed_file(fname, ht_exons)
        ga[fname] = _ga
    return ga, reads_by_gene, counts_by_gene


def read_a_bed_file(fname, ht_exons):
    import time
    start_time = time.time()
    counts = collections.defaultdict(int)
    reads_by_gene = collections.defaultdict(set)
    _ga = HTSeq.GenomicArray('auto', stranded=True)
    if not os.path.exists(fname):
        print "bed file not found: {o}".format(o=fname)
        return counts, reads_by_gene
    for n, line in enumerate(open(fname)):
        # if n > 1e6: break
        if not n % 1e6:
            min_per_line = float(time.time() - start_time)/float(60 * max([1., n]))
            print "Line {n}: Elapsed: {te} s. Minutes for 10 million reads: {f}".format(
                n=n, te=float(time.time() - start_time),
                f=min_per_line * 1e7/60.
            )
        s = line.rstrip('\n').split('\t')
        read_iv = HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), s[5])
        _ga[read_iv] += 1
        gene_ids = set()
        for iv, val in ht_exons[read_iv].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[gene_id] += 1
            reads_by_gene[gene_id].add(read_iv)
    print "Stats on {a}:".format(a=fname)
    cluster_reads.stats(reads_by_gene)
    return _ga, reads_by_gene, counts


def dir_creation(lib):
    dir_list = ['bedgraphs/', 'bed_uncollapsed/' 'bed_collapsed/',
                'figs/', 'logs/']
    for adir in dir_list:
        if not os.path.exists(lib['top'] + '/' + adir):
            os.system('mkdir ' + lib['top'] + '/' + adir)
        if not os.path.exists(lib['clusters_dir']):
            os.system('mkdir ' + lib['clusters_dir'])


if __name__ == '__main__':
    args = read_args()
    print args
    sys.path.insert(0, args.config_dir)
    import config
    lib = config.config()
    del sys.path[0]
    for key in lib:
        print key + ':' + lib[key]
    dir_creation(lib)
    logging.basicConfig(
        filename='logs/%s_find_peaks_by_perm.log' % datetime.datetime.now().strftime('%Hh%Mm'),
        level=logging.DEBUG)
    logging.info('Module called %s' % str(time.localtime()))
    logger = logging.getLogger(__name__)
    # Load raw reads.
    gtf = HTSeq.GFF_Reader(lib['gtf_raw'], end_included=True)
    reads_by_gene, counts_by_gene = get_bed_files(
        bed_folder=lib['read_beds'], gtf=gtf,
        args=args, lib=lib)
    # (sequences, rc_sequences, chr_lens) = subpeaks.get_sequences(lib)
    # (utr, exon, cds, other) = determine_feature_locations.get_gtf(args, lib=lib)
    # get_txpts() uses the peaks argument to select what txpts to load by a call to
    #  set(peaks['gene_name'].tolist()). That is, peaks is a dataframe.
    # peaks = pandas.DataFrame([{'gene_name': x} for x in cds.keys()])
    # print cds.keys()[:100]
    # (sequences, rc_sequences, chr_lens) = determine_feature_locations.get_sequences()
    # determine_feature_locations.flip_minus_strand_features(
    #         sequences, rc_sequences, chr_lens, peaks, utr, exon, cds, other, no_peaks=True)
    # txpts = determine_feature_locations.get_txpts(
    #     args, sequences, rc_sequences, peaks, utr, exon, cds, other, lib=lib
    # )
    exons_as_rows = p_values_of_clusters.get_exonic_ranges(lib['gtf_with_names'])
    while True:
        try:
            cluster_reads.run(
                reads_by_gene, counts_by_gene, exons_as_rows,
                lib=lib, args=args)
            print "Successfully finished."
        except:
            print traceback.format_exc()
            print "Failed."
        print "Hit enter to reload, CTRL-C to close."
        sys.stdin.readline()
        reloaded = False
        while not reloaded:
            try:
                reload(config)
                reload(cluster_reads)
                reload(p_values_of_clusters)
                reload(determine_feature_locations)
                reload(subpeaks)
                reload(put_reads_in_gene)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()

