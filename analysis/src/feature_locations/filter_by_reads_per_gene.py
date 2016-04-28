import cluster_combine

import pandas
import HTSeq
import argparse
import os
import sys
import numpy as np
import pandas
import collections
import glob
import re


def sum_reads_in_gene(full):
    for i, row in full.iterrows():
        full.loc[i, 'exp'] = row['fog_GGCA_counts.txt'] + row['fog_TGGC_counts.txt'] +\
            row['fog_GGTT_counts.txt'] + row['fog_CGGA_counts.txt']
        full.loc[i, 'control'] = row['control_CCGG_counts.txt'] + row['control_AATA_counts.txt'] +\
            row['control_TTGT_counts.txt'] + row['control_TTAA_counts.txt']
        full.loc[i, 'ratio'] = float(full.loc[i, 'exp'])/max([1., full.loc[i, 'control']])
    # for i, row in full.iterrows():
    #     full.loc[i, 'exp'] = row['fog_GGCA_reads'] + row['fog_TGGC_reads'] +\
    #         row['fog_GGTT_reads'] + row['fog_CGGA_reads']
    #     full.loc[i, 'control'] = row['control_CCGG_reads'] + row['control_AATA_reads'] +\
    #         row['control_TTGT_reads'] + row['control_TTAA_reads']
    #     full.loc[i, 'ratio'] = float(full.loc[i, 'exp'])/max([1., full.loc[i, 'control']])
    return full


def apply_cutoff(full, ratio_cutoff=10, min_reads_cutoff=200,
                 only_mrna=False):
    past = full[full['ratio']>=ratio_cutoff]
    past = past[past['exp']>min_reads_cutoff]
    if only_mrna:
        past = past[past['biotype']=='protein_coding']
    return past


if __name__ == '__main__':
    fname = sys.argv[1]
    table = pandas.read_csv(fname, sep='\t')
    table = cluster_combine.get_rpkm(table, {'gtf': 'lib/gtf_with_names_column.txt'})
    table = sum_reads_in_gene(table)
    print 'ratio: min {z} max {zz} mean {aa}'.format(
        z=min(table['ratio'].tolist()),
        zz=max(table['ratio'].tolist()),
        aa=np.mean(table['ratio'].tolist())
    )
    print 'exp: min {z} max {zz} mean {xx}'.format(
        z=min(table['exp'].tolist()),
        zz=max(table['exp'].tolist()),
        xx=np.mean(table['exp'].tolist()),
    )
    table = apply_cutoff(table)
    table.to_csv('tables/past_cutoff.txt', sep='\t')
    # print table.iloc[0]

