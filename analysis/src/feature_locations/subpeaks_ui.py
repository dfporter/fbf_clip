import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys

import subpeaks
import determine_feature_locations
import flocs
import build_image_objects_for_heatmaps
import peak_in_gene_line_plot
import scatterplot_correlation_by_wig
import output_heatmap_figures
import build_peaks_arr
import argparse
import pulsar
import edge_of_txpt
import config


def read_args():
    parser = argparse.ArgumentParser(description='''
        Determine the locations of features in genes across the genome.''')
    parser.add_argument(
        '-n', '--no_user_input', action='store_true', default=False,
        help='Do not request user input or run as a loop (Default: False).'
    )
    parser.add_argument('-i', '--input',
                        default=None,
                        help='''Input peaks file (Required).''')
    parser.add_argument('-c', '--config_ini',
        help='config.ini file',
        default='/groups/Kimble/Aman Prasad/clip/analysis/src/config.ini')
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
    lib = config.config(filepath=args.config_ini)
    # Logging.
    if not os.path.exists('logs'): os.system('mkdir logs')
    logging.basicConfig(
        filename='logs/%s_subpeaks.log' % datetime.datetime.now().strftime('%Hh%Mm'),
        level=logging.DEBUG)
    logging.info('Module called %s' % str(time.localtime()))
    logger = logging.getLogger(__name__)
    if not os.path.isfile(args.input):
        raise IOError('Input peaks file (-i) was not found.')
#    input_dir = os.path.dirname(args.input) + '/' #'../clip/pypeaks_fdr1_negip_local/'
    print "Loading bedgraphs..."
    ga = subpeaks.get_bedgraph(
        bedgraphs_folder=lib['coverage_wigs'],
        args=args, lib=lib)
    (peaks, txpts, chr_lens) = determine_feature_locations.get_data(args, lib=lib)
    (sequences, rc_sequences, chr_lens) = subpeaks.get_sequences(lib)
    print "Running script..."
    if args.no_ui:
        subpeaks.run(
            args.input, ga, lib, txpts, sequences, rc_sequences,
            chr_lens, args=args)
        sys.exit()
    while True:
        try:
            subpeaks.run(
                args.input, ga, lib, txpts, sequences, rc_sequences,
                chr_lens, args=args)
            print "Successfully finished."
        except:
            print traceback.format_exc()
            print "Failed."
        print "Hit enter to reload, CTRL-C to close."
        sys.stdin.readline()
        reloaded = False
        while not reloaded:
            try:
                reload(subpeaks)
                reload(config)
                reload(flocs)
                reload(build_image_objects_for_heatmaps)
                reload(peak_in_gene_line_plot)
                reload(scatterplot_correlation_by_wig)
                reload(output_heatmap_figures)
                reload(determine_feature_locations)
                reload(build_peaks_arr)
                reload(pulsar)
                reload(edge_of_txpt)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()
