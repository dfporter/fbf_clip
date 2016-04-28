__author__ = 'dp'
import pandas
import scipy
import scipy.stats
import numpy as np
import HTSeq
import sys
import re
import os


def load_bedgraph(bedgraphs_folder='data/wigs_coverage/'):
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    with open(bedgraphs_folder + 'both_fbfs_plus.bed', 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] = float(s[3])
    with open(bedgraphs_folder + 'both_fbfs_minus.bed', 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] = float(s[3])
    return ga


def get_max(ga, peak_iv):
    max_coverage = 0
    for iv, score in ga[peak_iv].steps():
        max_coverage = max(max_coverage, score)
    return max_coverage

def get_iv_from_pat(motif, row):
    _fbe_positions = []
    for fbe_pos in re.finditer(motif, row['seq']):
        if row['strand'] == '+':
            _fbe_positions.append(HTSeq.GenomicInterval(
                row['chrm'],
                row['left'] + fbe_pos.start(),
                row['left'] + fbe_pos.end(),
                row['strand']))
        else:
            flipped = [row['left'] + (len(row['seq']) - fbe_pos.end()),
                       row['left'] + (len(row['seq']) - fbe_pos.start())]
            _fbe_positions.append(HTSeq.GenomicInterval(
                row['chrm'], flipped[0], flipped[1], row['strand']
            ))
    return _fbe_positions


if __name__ == '__main__':
    df = pandas.read_csv(sys.argv[1], sep='\t')
    at_peaks = []
    at_fbes = []
    at_dhhs = []
    at_seven_mers = []
    at_minus_twos = []
    peak_positions = []
    fbe_positions = []
    dhh_positions = []
    seven_mers = []
    minus_twos = []
    ga = load_bedgraph()
    # Read motif ranges from peaks file
    for index, row in df.iterrows():
        peak_iv = HTSeq.GenomicInterval(
            row['chrm'], row['left'], row['right'], row['strand'])
        peak_positions.append(peak_iv)
        _fbes = get_iv_from_pat('tgt\w\w\wat', row)
        fbe_positions.extend(_fbes)
        minus_twos.extend(get_iv_from_pat('c\wtgt\w\w\wat', row))
        dhh_positions.extend(get_iv_from_pat('tgt[^c][^g][^g]at', row))
        _cryptics = get_iv_from_pat('tgt[ag][ac]at', row)
        seven_mers.extend(_cryptics)
    # Get the max signal in these ranges.
    for iv in peak_positions:
        at_peaks.append(get_max(ga, iv))
    for iv in minus_twos:
        at_minus_twos.append(get_max(ga, iv))
    for iv in fbe_positions:
        at_fbes.append(get_max(ga, iv))
    for iv in seven_mers:
        at_seven_mers.append(get_max(ga, iv))
    for iv in dhh_positions:
        at_dhhs.append(get_max(ga, iv))
    # Output stats.
    print "Median max coverage in a peak:"
    print np.median(np.array(at_peaks))
    print "Median max coverage at a -2C FBE (CNUGUNNNAU) in a peak:"
    print np.median(np.array(at_minus_twos))
    print "Median max coverage at an FBE (UGUNNNAU) in a peak:"
    print np.median(np.array(at_fbes))
    print "Median max coverage at an FBE (UGUDHHAU) in a peak:"
    print np.median(np.array(at_dhhs))
    print "Median max coverage at a 7-mer (UGU[AG][AC]AU) in a peak:"
    print np.median(np.array(at_seven_mers))
    print "t-test for fbe (UGUNNNAU) vs 7-mer:"
    print scipy.stats.ttest_ind(np.array(at_fbes), np.array(at_seven_mers))
    print "t-test for fbe (UGUDHHAU) vs 7-mer:"
    print scipy.stats.ttest_ind(np.array(at_dhhs), np.array(at_seven_mers))
    print "t-test for UGUDHHAU vs UGUNNNAU:"
    print scipy.stats.ttest_ind(np.array(at_dhhs), np.array(at_fbes))
