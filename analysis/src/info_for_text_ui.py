import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys
import argparse

import info_for_text
import add_height
import compare_with_ripchip
from load_bedgraph import load_bedgraph


def get_bedgraph():
    ga = {}
    #ga['both'] = load_bedgraph('extra_bedgraph_unnorm/combined_fbf.txt')
    ga['fbf1'] = load_bedgraph('bedgraph_unnorm/combined_fbf1')
    ga['both'] = load_bedgraph('bedgraph_unnorm/combined_fbf')
    ga['fbf2'] = load_bedgraph('bedgraph_unnorm/combined_fbf2')
    return ga


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--five_reps',
        default='/groups/Kimble/Aman Prasad/clip/methods/filtered_gauss_nil_01_five_reps/',
        help='''Input directory of peaks, 5/6 reps.''')
    parser.add_argument(
        '-s', '--separate',
        default='/groups/Kimble/Aman Prasad/clip/methods/filtered_gauss_same_n2_nil_01/',
        help='''Input directory of peaks, separate.''')
    args = parser.parse_args()
    # Logging.
    ga = get_bedgraph()
    print "Running script..."
    while True:
        try:
            info_for_text.add_info_for_text(ga, args)
            print "Successfully finished."
        except:
            print traceback.format_exc()
            print "Failed."
        print "Hit enter to reload, CTRL-C to close."
        sys.stdin.readline()
        reloaded = False
        while not reloaded:
            try:
                reload(add_height)
                reload(info_for_text)
                reload(compare_with_ripchip)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()
