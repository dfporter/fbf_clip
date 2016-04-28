__author__ = 'dp'

import pandas
import re
import operator
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle
import matplotlib
import copy
import build_peaks_arr
import scipy
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import linkage

def by_location(peaks, output_dirname):
    cds = peaks[peaks['location']=='CDS']
    fiveprime = peaks[peaks['location']=="""5'UTR"""]
    threeprime = peaks[peaks['location']=="""3'UTR"""]
    with open(output_dirname + 'cds_peaks.fa','w') as f:
        n = 0
        for index, row in cds.iterrows():
            n += 1
            f.write('>{id}_{n}\n{seq}\n'.format(
                id=cds.loc[index, 'gene_name'], n=n, seq=cds.loc[index, 'seq']))
    with open(output_dirname + 'fiveprime_peaks.fa', 'w') as f:
        n = 0
        for index, row in fiveprime.iterrows():
            n += 1
            f.write('>{id}_{n}\n{seq}\n'.format(
                id=fiveprime.loc[index, 'gene_name'], n=n, seq=fiveprime.loc[index, 'seq']))
    with open(output_dirname + 'threeprime_peaks.fa', 'w') as f:
        n = 0
        for index, row in threeprime.iterrows():
            n += 1
            f.write('>{id}_{n}\n{seq}\n'.format(
                id=threeprime.loc[index, 'gene_name'], n=n, seq=threeprime.loc[index, 'seq']))

def highest_vs_secondary(peaks, output_dirname):
    highest_peak_seqs = []
    secondary_peak_seqs = []
    j = 0
    for gene_name in set(peaks['gene_name'].tolist()):
        j += 1
        n = 0
        peaks_in_gene = peaks[peaks['gene_name']==gene_name]
        if len(peaks_in_gene.index) == 1:
            highest_peak_seqs.append(
                ('gene_%i_peak_%i' % (j, n), peaks_in_gene.iloc[0].seq))
        elif len(peaks_in_gene.index) > 1:
            seqs = peaks_in_gene['seq'].tolist()
            heights = peaks_in_gene['height'].tolist()
            by_height = sorted(zip(seqs, heights), key=lambda x: x[1], reverse=True)
            highest_peak_seqs.append(
                ('gene_%i_peak_%i' % (j, n), peaks_in_gene.iloc[0].seq)
            )
            #for index in range(1, len(peaks_in_gene.index)):
            for i, tup in enumerate(by_height[1:]):
                secondary_peak_seqs.append(
                    ('gene_%i_peak_%i' % (j, i), tup[0])
                )
    with open(output_dirname +'/highest_peak_seqs.fa', 'w') as f:
        for tup in highest_peak_seqs:
            #id = re.sub('\s', '_', str(tup[0])).partition('Name')[0]
            f.write(">{id}\n{seq}\n".format(id=tup[0], seq=tup[1]))
    with open(output_dirname +'/secondary_peak_seqs.fa', 'w') as f:
        for tup in secondary_peak_seqs:
            #id = re.sub('\s', '_', str(tup[0])).partition('Name')[0]
            f.write(">{id}\n{seq}\n".format(id=tup[0], seq=tup[1]))

if __name__ =='__main__':
    parser = argparse.ArgumentParser(description='''
        Get peak sequences, as split by order of peak height or location in the gene.''')
    parser.add_argument('-i', '--input',
                        default='new_pypeaks_fdr1_negip_local_five_reps/combined_fbf.txt',
                        help='''Input peaks file (Required).''')
    args = parser.parse_args()
    args.input = args.input
    peaks = pandas.read_csv(args.input, sep='\t')
    output_dirname = 'data/%s/fasta/' % os.path.basename(os.path.dirname(args.input))
    print output_dirname
    if not os.path.exists('data/%s' % os.path.basename(os.path.dirname(args.input))):
        os.system('mkdir data/%s' % os.path.basename(os.path.dirname(args.input)))
    if not os.path.exists(output_dirname):
        os.system('mkdir ' + output_dirname)
    #sys.exit()
    highest_vs_secondary(peaks, output_dirname)
    by_location(peaks, output_dirname)