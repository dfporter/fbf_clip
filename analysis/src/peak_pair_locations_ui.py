import traceback
import re
import sys
import pandas
import numpy as np
import logging
import time
import os
import datetime
import argparse
import glob

def read_args():
    parser = argparse.ArgumentParser(description='''
    Call peaks in CLIP-seq data.
    Outputs first pickled objects to data/raw_peak_objs_by_chrm_from_callpeaks_*,
    and equivalent dataframe-based tables to data/cwt_calls/.
    ''')
    parser.add_argument(
        '-n', '--no_user_input', action='store_true', default=False,
        help='Do not request user input or run as a loop (Default: False).'
    )
    parser.add_argument(
        '-p', '--peaks', default='gauss_nil/combined_fbf.txt',
        help='Peaks filename (Required).'
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # Logging.
    if not os.path.exists('logs'): os.system('mkdir logs')
    logging.basicConfig(
        filename='logs/%s_peak_pair_locations.log' % datetime.datetime.now().strftime('%Hh%Mm'),
        level=logging.DEBUG)
    logging.info('Module called %s' % str(time.localtime()))
    logger = logging.getLogger(__name__)
    # Get arguments and run.
    args = read_args()
    import peak_pair_locations
    gtf_filename = 'lib/gtf_with_names_column.txt'
    gtf_df = pandas.read_csv(gtf_filename, sep='\t')
    print gtf_df.head()
    gtf_d = peak_pair_locations.convert_df_to_dict(gtf_df)
    if args.no_user_input:
        if os.path.isdir(args.peaks):
            for fname in glob.glob(args.peaks + '/*txt'):
                print "Acting on %s" % fname
                peak_pair_locations.run(fname, gtf_d)
        elif os.path.isfile(args.peaks):
            peak_pair_locations.run(args.peaks, gtf_d)
        else:
            print "Input %s was neither a file nor directory" % args.peaks
            raise Exception
        sys.exit()
    while True:
        try:
            if os.path.isdir(args.peaks):
                for fname in glob.glob(args.peaks + '/*txt'):
                    print "Acting on %s" % fname
                    peak_pair_locations.run(fname, gtf_d)
            elif os.path.isfile(args.peaks):
                peak_pair_locations.run(args.peaks, gtf_d)
            else:
                print "Input %s was neither a file nor directory" % args.peaks
                raise Exception
            print "Successfully finished."
        except:
            print traceback.format_exc()
            print "Failed."
        print "Hit enter to reload, CTRL-C to close."
        sys.stdin.readline()
        reloaded = False
        while not reloaded:
            try:
                reload(peak_pair_locations)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()
