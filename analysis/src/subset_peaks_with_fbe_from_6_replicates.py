import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import pysam
import matplotlib
import sys
from matplotlib_venn import *
import os
import HTSeq
from subset_peaks_with_fbe import *


def combine_peaks(replicates):
    combined = {}  # Key: (chrm, left, strand) tuple.
    first_rep_name = replicates.keys()[0]
    as_single_df = replicates[first_rep_name]
    as_single_df = pandas.concat([replicates[rep][:] for rep in replicates], ignore_index=True)
    print "after append: size of combined df: %i." % len(as_single_df)
    for index, row in as_single_df.iterrows():
        if index % 1000 == 0: print "combining %i/%i." % (index, len(as_single_df))
        gene = as_single_df.loc[index, 'gene_name']
        iv = [as_single_df.loc[index, 'chrm'], as_single_df.loc[index, 'left'],
              as_single_df.loc[index, 'right'], as_single_df.loc[index, 'strand']]
        overlapping_peaks = overlapping(iv, as_single_df)
        if len(overlapping_peaks) > 3:
            consensus = consensus_peak(overlapping_peaks)
            tup = (consensus['chrm'], consensus['left'], consensus['strand'])
            if tup not in combined:
                combined[tup] = consensus
    return combined


def consensus_peak(peaks):
    rows = peaks[peaks['height']==max(peaks['height'])]
    return rows.iloc[0].to_dict()


def overlapping(iv, comb_df):
    overlap = comb_df[(comb_df['chrm']==iv[0]) & (comb_df['strand']==iv[3]) & (comb_df['left']<=iv[2]) & (iv[1]<=comb_df['right'])]
    return overlap


if __name__ == '__main__':
    top_level_dir = sys.argv[1]
    replicates = {}
    for filename in glob.glob(top_level_dir + '/fbf*'):
        replicates[filename] = pandas.read_csv(filename, sep='\t')
    combined_dict = combine_peaks(replicates)
    combined_df = pandas.DataFrame(combined.values())
    combined_df.to_csv('with_fbe_combined_fbfs_%s/%s' % (
        top_level_dir, 'fbf1_and_2'), sep='\t')

