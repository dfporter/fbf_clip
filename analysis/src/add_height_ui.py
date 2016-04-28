import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys

def get_bedgraph(
        do_combine_bedgraphs=False, bedgraphs_folder='data/wigs_five_prime/',
        lib=None):
    plus_file = bedgraphs_folder + 'both_fbfs_plus.bed'
    minus_file = bedgraphs_folder + 'both_fbfs_minus.bed'
    ga = {}
    ga['both'] = get_a_bedgraph(plus_file, minus_file)
    plus_file = bedgraphs_folder + 'fbf1_reads_plus.bed'
    minus_file = bedgraphs_folder + 'fbf1_reads_minus.bed'
    ga['fbf1'] = get_a_bedgraph(plus_file, minus_file)
    plus_file = bedgraphs_folder + 'fbf2_reads_plus.bed'
    minus_file = bedgraphs_folder + 'fbf2_reads_minus.bed'
    ga['fbf2'] = get_a_bedgraph(plus_file, minus_file)
    return ga

def get_a_bedgraph(plus_file, minus_file):
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    with open(plus_file, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] = float(s[3])
    with open(minus_file, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] = float(s[3])
    return ga
import add_height

if __name__ == '__main__':
    # Logging.
    ga = get_bedgraph(
        bedgraphs_folder='/groups/Kimble/Aman Prasad/clip/data/wigs_coverage/')
    print "Running script..."
    while True:
        try:
            add_height.add_height(ga)
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
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()
