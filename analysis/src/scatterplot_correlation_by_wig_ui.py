import sys
import glob
import os
import re
import pandas
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
import pickle
import argparse

sys.path.insert(0, '/Network/Servers/file.biochem.wisc.edu/Volumes/BioSAN/Users/dfporter/.local/lib/python2.6/site-packages')
# Need HTSeq version 0.6, hence the sys.path alteration.


import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys

import add_height
import scatterplot_correlation_by_wig
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Add sequence, has_fbe, rip_column, ect.
        info columns to a given directory -d of peaks files.
        Can also combine replicates and output the combined list.
        Will apply a normalization ratio of ~1.6 to FBF-2 reads.'''
    )
    parser.add_argument('-d', '--directory', dest='directory',
                        default='new_pypeaks_fdr1_negip_local/',
                        help='Input directory of peaks files.')
    parser.add_argument('-p', '--load_data_from_pickle',
                        default=False,
                        action='store_true',
                        help='''Load data from the file \
                        data/integrals_from_scatterplot_correlation_by_wig.p''')
    parser.add_argument('-n', '--no_ui',
                        default=False,
                        action='store_true',
                        help='''Don't expect any user input (default: False).''')
#    parser.add_argument('-c', '--cols', default=False, action='store_true',
#                        help='''Use the reads_fbf1/reads_fbf2 columns from the \
#                        input peaks file, as generated by \
#                        filter_peaks_by_ratio_in_peak_region.py --load_bed_info''')
    args = parser.parse_args()
    ga = scatterplot_correlation_by_wig.get_ga_from_bedgraph(args)
    print "Running script..."
    if args.no_ui:
        scatterplot_correlation_by_wig.run(args, ga)
        sys.exit()
    while True:
        try:
            scatterplot_correlation_by_wig.run(args, ga)
            print "Successfully finished."
        except:
            print traceback.format_exc()
            print "Failed."
        print "Hit enter to reload, CTRL-C to close."
        sys.stdin.readline()
        reloaded = False
        while not reloaded:
            try:
                reload(scatterplot_correlation_by_wig)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()
