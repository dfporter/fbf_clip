"""
Examines the actual ratios between FBF-1 and FBF-2 in cases
in which the peak caller has called a gene uniquely for one or the other.

Takes input of a /counts directory with reads/gene counts for -1 or -2.
Takes input of a peaks list for -1 or -2.
gtf[raw]

Outputs a set of overlap numbers at different min_read cutoffs.
Outputs a histogram of read/gene ratios.

"""
__author__ = 'dp'

import os
import sys
import re
import glob

import pandas
import matplotlib.pyplot as plt


def read_as_table(fname, use_col='WB ID', use_header_val=None):
    _table = {}
    key = use_header_val
    df = pandas.read_csv(fname, sep='\t')
    dfd = df.to_dict('records')
    for row in dfd:
        _table[row[key]] = row
    return _table


def add_wb_name(
        counts, gtf, fname='./lib/gtf_with_id_column.txt'):
    print "Adding wb name..."
    print str(gtf)[:100]
    # countsl = counts.to_dict('records')
    # countsd = {(row['gene_name'], row)}
    counts['gene_id'] = counts['gene_name']
    for i, row in counts.iterrows():
        if row['gene_name'] in gtf:
            gene_name = gtf[row['gene_name']]['gene_name']
            counts.loc[i, 'gene_name'] = gene_name
        else:
            print 'gene name not in wb gtf: %s' % row['gene_name']
    return counts
    # cols = get_cols(peaks_list)
    # peaks_df.to_csv(peaks_fname, cols=cols, sep='\t', index=False)


def counts_as_dict(counts_s):
    rows = counts_s.to_dict('records')
    # print rows
    counts_sd = {}
    for row in rows:
        counts_sd[row['gene_name']] = row['reads']
    return counts_sd


def add_counts(peaks, counts_sd):
    peaks['reads_in_gene'] = 0
    for i, row in peaks.iterrows():
        if row['gene_name'] not in counts_sd:
            continue
        peaks.loc[i, 'reads_in_gene'] = counts_sd[row['gene_name']]
    return peaks


def get_overlaps_by_min_read_number(
        counts_s1, counts_s2, peaks1, peaks2, cutoff=10):
    """
    counts_s1/s2: dataframe
    counts_s1d/s2d: dict
    peaks1s/2s: dataframe
    """
    print "Adding counts..."
    counts_s1d = counts_as_dict(counts_s1)
    counts_s2d = counts_as_dict(counts_s2)
    peaks1 = add_counts(peaks1, counts_s1d)
    peaks2 = add_counts(peaks2, counts_s2d)
    print "Filtering..."
    peaks1s, peaks2s = filter_by(counts_s1, counts_s2, peaks1, peaks2,
                                 cutoff=cutoff)
    print "Getting ratios and making figures..."
    in_one, in_two, in_both, li = divide_into_categories(peaks1s, peaks2s)
    all_targs = in_one | in_two | in_both
    ratios_df = get_ratios(counts_s1d, counts_s2d, all_targs, return_df=True)
    ratios_df.to_csv('enrichment_ratios.txt', sep='\t')
    ratios_df.sort(columns=['ratio_fbf1_over_fbf2'], inplace=True)
    ratios = get_ratios(counts_s1d, counts_s2d, in_both)
    top_fbf2 = ratios_df.head(100)
    top_fbf2.to_csv('fbf2_enriched.txt', columns=['gene_name'], sep='\t', index=False)
    top_fbf1 = ratios_df.tail(100)
    top_fbf1.to_csv('fbf1_enriched.txt', columns=['gene_name'], sep='\t', index=False)
    plt.clf()
    f, ax = plt.subplots(2, 2)
    print li
    ax[0][0] = make_hist(ratios, ax[0][0], title='Joint targets', output_filename='joint_targets.pdf')
    ratios = get_ratios(counts_s1d, counts_s2d, in_one)
    print li
    ax[0][1] = make_hist(ratios, ax[0][1], title='FBF-1 only targets', output_filename='only_fbf1_targets.pdf')
    ratios = get_ratios(counts_s1d, counts_s2d, in_two)
    print li
    ax[1][0] = make_hist(ratios, ax[1][0], title='FBF-2 only targets', output_filename='only_fbf2_targets.pdf')
    output_filename = 'test.pdf'
    plt.tight_layout()
    plt.savefig(output_filename)
    plt.clf()
    plt.close()


def make_hist(ratios, ax, title='', output_filename='test.pdf'):
    ax.hist(ratios, 100, facecolor='black')
    ax.set_xlabel('log2(FBF-1/FBF-2) reads per gene')
    ax.set_ylabel('Number of RNAs')
    ax.set_title(title)
    # plt.savefig(
    #    output_filename, format='pdf')
    return ax


def write_outliers():
    pass


def get_ratios(counts1, counts2, gene_names, return_df=False):
    ratios = []
    ratios_d = []
    import numpy as np
    for name in gene_names:
        if name in counts1:
            numerator = float(counts1[name])
        else:
            numerator = 0.
        if name in counts2:
            denominator = float(counts2[name])
        else:
            denominator = 1.
        ratio = np.log2(numerator/denominator)
        ratios.append(ratio)
        ratios_d.append({'gene_name': name, 'ratio_fbf1_over_fbf2': ratio})
    ratios_d = pandas.DataFrame(ratios_d)
    ratios_d.sort(columns=['ratio_fbf1_over_fbf2'], inplace=True)
    print """
    Minimum ratios: {n}
    Return ratios: {a}
    """.format(
        n=ratios_d.head(),
        a=ratios_d.tail()
    )
    if return_df: return ratios_d
    return ratios


def filter_by(counts1, counts2, peaks1s, peaks2s, cutoff=10):
    counts_s1 = counts1[counts1['reads']>=cutoff]
    counts_s2 = counts2[counts2['reads']>=cutoff]
    above_cutoff = set(counts_s1['gene_name']) & set(counts_s2['gene_name'])
    peaks1s['above_cutoff'] = 0
    peaks2s['above_cutoff'] = 0
    for index, row in peaks1s.iterrows():
        # print 'is {a} in {b}'.format(a=row['gene_name'],
        #                              b=str(above_cutoff)[:100])
        if row['gene_name'] in above_cutoff:
            peaks1s.loc[index, 'above_cutoff'] = 1
    for index, row in peaks2s.iterrows():
        if row['gene_name'] in above_cutoff: peaks2s.loc[index, 'above_cutoff'] = 1
    _peaks1s = peaks1[peaks1s['above_cutoff']>0]
    _peaks2s = peaks2[peaks2s['above_cutoff']>0]
    return _peaks1s, _peaks2s


def divide_into_categories(peaks1s, peaks2s):
    in_one = set(peaks1s['gene_name'].tolist()) - set(peaks2s['gene_name'])
    in_two = set(peaks2s['gene_name'].tolist()) - set(peaks1s['gene_name'])
    in_both = set(peaks1s['gene_name'].tolist()) & set(peaks2s['gene_name'])
    li = """FBF-1 only: {n} FBF-2 only {a} In both: {v}""".format(
        n=len(in_one), a=len(in_two), v=len(in_both)
    )
    return in_one, in_two, in_both, li


if __name__ == '__main__':
    peaks_dir = sys.argv[1]
    counts_dir = sys.argv[2]
    peaks_file1 = peaks_dir + '/combined_fbf1.txt'
    peaks_file2 = peaks_dir + '/combined_fbf2.txt'
    counts_file1 = counts_dir + '/combined_fbf1_counts.txt'
    counts_file2 = counts_dir + '/combined_fbf2_counts.txt'
    peaks1 = pandas.read_csv(peaks_file1, sep='\t')
    peaks2 = pandas.read_csv(peaks_file2, sep='\t')
    peaks1_onlymrna = peaks1[peaks1['biotype']=='protein_coding']
    peaks2_onlymrna = peaks2[peaks2['biotype']=='protein_coding']
    li = """
    Peaks file 1: {a}, genes: {b}
    Peaks file 2: {c}, gene: {d}
    Genes in either target list: {z}
    Genes in either target list if you don't count ncRNA: {pc_only}
    Overlapping targets: {z2}
    Is ncRNA present in peaks file 1? {yn}""".format(
        a=len(peaks1.index), b=len(set(peaks1['gene_name'].tolist())),
        c=len(peaks2.index), d=len(set(peaks2['gene_name'].tolist())),
        z=len(set(peaks1['gene_name'].tolist()) | set(peaks2['gene_name'].tolist())),
        pc_only=len(set(peaks1_onlymrna['gene_name'].tolist()) | set(peaks2_onlymrna['gene_name'].tolist())),
        z2=len(set(peaks1['gene_name'].tolist()) & set(peaks2['gene_name'].tolist())),
        yn='no' if len(peaks1[peaks1['biotype']=='protein_coding'].index) == len(peaks1.index) else 'yes',
    )
    print peaks1['biotype'].value_counts()
    print peaks2['biotype'].value_counts()
    print li
    # Counts.
    counts1 = pandas.read_csv(counts_file1, sep='\t', header=None, names=['gene_name', 'reads'])
    counts2 = pandas.read_csv(counts_file2, sep='\t', header=False, names=['gene_name', 'reads'])
    total_1 = float(counts1['reads'].sum())
    print total_1
    total_2 = float(counts2['reads'].sum())
    print total_2
    print '()()()'
    counts1['reads'] = [1e6 * float(x)/total_1 for x in counts1['reads']]
    counts2['reads'] = [1e6 * float(x)/total_2 for x in counts2['reads']]
    print counts1['reads'].sum()
    print counts2['reads'].sum()
    print "Loading gtf..."
    fname = '/groups/Kimble/Common/fog_iCLIP/calls/lib/gtf_with_names_column.txt'
    gtf = read_as_table(fname, use_header_val='gene_id')
    counts1 = add_wb_name(counts1, gtf)
    counts2 = add_wb_name(counts2, gtf, fname='/groups/Kimble/Common/fog_iCLIP/calls/lib/gtf_with_names_column.txt')
    get_overlaps_by_min_read_number(
        counts1, counts2, peaks1, peaks2, cutoff=10)



