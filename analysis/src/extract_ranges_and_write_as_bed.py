import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import HTSeq
import pysam
import sys
import matplotlib
import os
import xlwt


class peak:
    def __init__(self, chrm, left, right, strand):
        self.chrm = str(chrm)
        self.left = int(left)
        self.right = int(right)
        self.strand = str(strand)
        self.iv = HTSeq.GenomicInterval(chrm, left, right, strand)
        self.local = {}
        self.exons = {}
        self.pvalues = {}

def merge_overlapping_on_chrm_and_strand(intervals):
    """Merge in a given chrom and strand.
    """
    sorted_by_lower_bound = sorted(intervals, key=lambda x: x.left)
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher.left <= lower.right:
                upper_bound = int(max(lower.right, higher.right))
                new_peak = peak(lower.chrm, lower.left, upper_bound, lower.strand)
                merged[-1] = new_peak  # replace by merged interval
            else:
                merged.append(higher)
    return merged

def merge_overlapping(joined_peaks):
    intervals = []
    for index, row in joined_peaks.iterrows():
        intervals.append(peak(
            joined_peaks.loc[index, 'chrm'],
            joined_peaks.loc[index, 'left'],
            joined_peaks.loc[index, 'right'],
            joined_peaks.loc[index, 'strand']))
    by_chrm = {}
    merged = {}
    for apeak in intervals:
        by_chrm.setdefault(apeak.chrm, {})
        merged.setdefault(apeak.chrm, {})
        by_chrm[apeak.chrm].setdefault(apeak.strand, [])
        merged[apeak.chrm].setdefault(apeak.strand, [])
        by_chrm[apeak.chrm][apeak.strand].append(apeak)
    for chrm in by_chrm:
        for strand in by_chrm[chrm]:
            merged[chrm][strand] = merge_overlapping_on_chrm_and_strand(by_chrm[chrm][strand])
            # Check.
            check_overlap(merged[chrm][strand])
    return merged


def check_overlap(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda x: x.left)
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher.left <= lower.right:
                print "Overlapping peak ranges..\n%s\n%s\n*" % (str(lower), str(higher))



def write_peak_ranges_as_bed(peaks, bed_filename):
    bed_file = open(bed_filename, 'a')
    for chrm in peaks:
        for strand in peaks[chrm]:
            for apeak in peaks[chrm][strand]:
                bed_file.write('{chrm}\t{left}\t{right}\t.\t1\n'.format(
                    chrm=apeak.chrm,
                    left=apeak.left,
                    right=apeak.right))
    bed_file.close()


if __name__ == '__main__':
    #top_level_dir = './pypeaks_fdr1_negip_local/'
    top_level_dir = sys.argv[1]
    peaks = {}
    # Clear the output file.
    bed_filename = './data/%s.bed' % os.path.basename(top_level_dir.rstrip('/'))
    if not os.path.exists('./data'): os.system('mkdir data')
    if not os.path.exists(bed_filename):
        os.system('mkdir ' + os.path.dirname(bed_filename))
    os.system('rm %s' % bed_filename)
    for filename in glob.glob(top_level_dir + '/combined_*.txt'):
        print filename
        os.system('wc -l %s' % filename)
        peaks[filename] = pandas.read_csv(filename, sep='\t')
    joined_peaks = pandas.concat([peaks[exp][:] for exp in peaks], ignore_index=True)
    amerged = merge_overlapping(joined_peaks)
    write_peak_ranges_as_bed(amerged, bed_filename)

