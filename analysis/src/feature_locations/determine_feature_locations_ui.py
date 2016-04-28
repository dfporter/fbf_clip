import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys
import config

import determine_feature_locations
import flocs
import build_image_objects_for_heatmaps
import peak_in_gene_line_plot
import scatterplot_correlation_by_wig
import output_heatmap_figures
import build_peaks_arr

def read_args():
    parser = argparse.ArgumentParser(description='''
        Determine the locations of features in genes across the genome.''')
    parser.add_argument(
        '-n', '--no_user_input', action='store_true', default=False,
        help='Do not request user input or run as a loop (Default: False).'
    )
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
    parser.add_argument('-c', '--config_ini',
        help='config.ini file',
        default='/groups/Kimble/Aman Prasad/clip/analysis/src/config.ini')
#    parser.add_argument('-w', '--continuous_signal',
#                        default=False, action='store_true',
#                        help='Load true mapping data to display peak signal (slow).')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = read_args()
    lib = config.config(filepath=args.config_ini)
    print lib
    # Logging.
    if not os.path.exists('logs'): os.system('mkdir logs')
    logging.basicConfig(
        filename='logs/%s_determine_feature_locations.log' % datetime.datetime.now().strftime('%Hh%Mm'),
        level=logging.DEBUG)
    logging.info('Module called %s' % str(time.localtime()))
    logger = logging.getLogger(__name__)
    (peaks, txpts, chr_lens) = determine_feature_locations.get_data(args, lib=lib)
    input_dir = os.path.dirname(args.input) + '/' #'../clip/pypeaks_fdr1_negip_local/'
    print "Running script..."
    while True:
        try:
            determine_feature_locations.make_figs(peaks, txpts, chr_lens, args, input_dir)
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
                reload(flocs)
                reload(build_image_objects_for_heatmaps)
                reload(peak_in_gene_line_plot)
                reload(scatterplot_correlation_by_wig)
                reload(output_heatmap_figures)
                reload(determine_feature_locations)
                reload(build_peaks_arr)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()
