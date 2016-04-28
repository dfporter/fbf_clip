"""
Load the clusters output by find_peaks_by_permutations.

"""
import HTSeq
import argparse
import traceback
import sys

import determine_feature_locations
import subpeaks
import put_reads_in_gene
import cluster_reads
import p_values_of_clusters
import cluster_combine


def read_args():
    parser = argparse.ArgumentParser(description='''
        Load clusters output by find_peaks_by_permutations.py
        and identify clusters overlapping (at significance) in multiple
        replicates.''')
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



if __name__ == '__main__':
    args = read_args()
    print args
    sys.path.insert(0, args.config_dir)
    import config
    lib = config.config()
    del sys.path[0]
    for key in lib:
        print key + ':' + lib[key]
    # Load raw reads.
    gtf = HTSeq.GFF_Reader(lib['gtf_raw'], end_included=True)
    # get_txpts() uses the peaks argument to select what txpts to load by a call to
    #  set(peaks['gene_name'].tolist()). That is, peaks is a dataframe.
    exons_as_rows = p_values_of_clusters.get_exonic_ranges(lib['gtf_with_names'])
    while True:
        try:
            cluster_combine.run(args, lib, gtf, exons_as_rows)
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
                reload(cluster_combine)
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

