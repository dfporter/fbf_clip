import HTSeq
import logging
# import datetime
import argparse
import os
# import time
import traceback
import sys
import numpy as np
import pandas
import collections
import glob
import re
import copy

import cluster_reads
import p_values_of_clusters


def get_ivs_by_gene(input_files):
    ivs_by_gene = collections.defaultdict(dict)
    for fname in input_files:
        print fname
        db = pandas.read_csv(fname, sep='\t')
        db = db[db['padj'] < 0.01]
        db = db[db['max_coverage']>1]
        if ('gene_name' not in db.columns) and ('gene_id' in db.columns):
            db['gene_name'] = db['gene_id']
        db = db.groupby('gene_name')
        for gene_name, rows in db:
            pvalues = rows['padj'].tolist()
            rows_list = [(x['chrom'], x['left'], x['right'], x['strand']) for i,x in rows.iterrows()]
            ivs_by_gene[gene_name][fname] = [HTSeq.GenomicInterval(*x) for x in rows_list]
    return ivs_by_gene


def get_reproducible_peaks(ivs_by_gene):
    clusters_by_gene = collections.defaultdict(list)
    for gene in ivs_by_gene:
        rig = []
        for fname in ivs_by_gene[gene]:
            rig.extend(ivs_by_gene[gene][fname])
        clusters = p_values_of_clusters.find_clusters_in_gene(rig)
        clusters = [x + [0] for x in clusters]
        for row in clusters:
            if row[2] > 4:  # Number of replicates with cluster.
                print clusters
        # clusters is a list of the following for each cluster in this gene:
        # [left, right,, max_coverage=number of replicates with peak,
        #       a zero for p value]
        filtered_clust = [clust for clust in clusters if clust[2] > 1]
        # for clust in clusters:
        #     if clust[2] > 2:
        #         filtered_clust.append(clust)
        # clusters = [x for x in clusters if x[2]>2]
        if len(filtered_clust) > 0:
            clusters_by_gene[gene] = filtered_clust #[x.extend([0]) for x in clusters]
    return clusters_by_gene


def scribe(cmd):
    print cmd
    os.system(cmd)


def get_combined_clusters_from_file_list(file_list, exons_as_rows, output_filename='temp.txt'):
    ivs_by_gene = get_ivs_by_gene(file_list)
    clusters_by_gene = get_reproducible_peaks(ivs_by_gene)
    # print clusters_by_gene
    if len(clusters_by_gene) > 0:
        table = cluster_reads.convert_to_table(clusters_by_gene, exons_as_rows)
        # print table['gene_id'].value_counts()
        # print len(table['gene_id'].value_counts())
        if not os.path.exists(os.path.dirname(output_filename)):
            scribe('mkdir ' + os.path.dirname(output_filename))
        try:
            table = get_counts_per_gene(table)
        except:
            print "TO DO: ADD GET_COUNTS_PER_GENE FUNCTION TO CLUSTER_COMBINE"
        cluster_reads.write_table(table, output_filename=output_filename)
        return table
    return pandas.DataFrame()


def verification(db):
    if 'max_coverage' not in db.columns:
        print "cluster_combine.verification(): No max coverage column in dataframe:"
        print db.head(1)
        return
    too_low = db[db['max_coverage']<2]
    print 'max coverage < 2: {i}'.format(i=len(too_low.index))
    print 'Peaks: {p}. Genes: {g}.'.format(
        p=len(db.index), g=len(set(db['gene_id']))
    )


def get_file_list(lib):
    file_list = [lib[x] for x in lib if (
        re.match('exp_bed.*', x) or (re.match('control_bed.*', x))
    )]
    return file_list


def get_gene_len(genes, filename, use='gene_id'):
    gtf = read_as_table_of_lists(filename, use_header_val='gene_id')
#    with open(filename, 'r') as f:
#        for li in f:
#            s = li.rstrip('\n').split('\t')
    lens = {}
    for gene in genes:
        if gene not in gtf:
            lens[gene] = 1e3
            print 'Missing ' + gene
            continue
        txpts = set([x['transcript_id'] for x in gtf[gene]]) - set([np.nan])
        try:
            this_txpt = list(txpts)[0]  # Expect only one txpt, but just in case.
        except:
            print gene
            print txpts
            print gtf[gene]
        lens[gene] = sum([
            int(row['4']) - int(row['3']) for row in gtf[gene] if (
                (row['transcript_id'] == this_txpt) and (row['2'] == 'exon')
            )
        ])
        if lens[gene] == 0:
            print this_txpt
            print set([x['transcript_id'] for x in gtf[gene]])
            print [x['2'] for x in gtf[gene]]
    return lens


def read_as_table_of_lists(fname, use_col='WB ID', use_header_val=None):
    if use_header_val is not None:
        _table = collections.defaultdict(list)
        key = use_header_val
        df = pandas.read_csv(fname, sep='\t')
        dfd = df.to_dict('records')
        for row in dfd:
            _table[row[key]].append(row)
        return _table


def get_rpkm(
        table, lib,
        counts_fname='counts/combined_counts.txt'):
    db = pandas.read_csv(counts_fname, sep='\t', index_col=False)
    if ('gene' in db.columns) and ('gene_id' not in db.columns):
        db['gene_id'] = db['gene']
    lens = get_gene_len(db['gene'].tolist(), lib['gtf'], use='gene_id')
    lens_as_list = [float(lens[x]) for x in db['gene'].tolist()]
    read_cols = [x for x in db.columns if x[:4] != 'gene']
    for col in read_cols:
        total = float(db[col].sum())
        db[col] = [1e6 * float(x)/total for x in db[col].tolist()]
        db[col] = [1e3 * float(x)/max([1, lens_as_list[i]]) for i, x in enumerate(db[col].tolist())]
    print """cluster_combine.get_rpkm():
read_cols from file {a}: {s}. Top line after conversion from raw reads to rpkm: {k}
***""".format(
        a=counts_fname, s=read_cols, k=db.iloc[0]
    )
    full = pandas.merge(db, table, on='gene_id')
    return full


def run(args, lib, gtf, exons_as_rows):
    bed_files = get_file_list(lib)
    cluster_file_list = []
    for bed_file in bed_files:
        cluster_file_list.append(
            lib['clusters_dir'] + '/' + os.path.basename(
                bed_file.rstrip('.bed')) + '.txt'
        )
    print "cluster_combine.run() combining clusters from files:"
    print cluster_file_list
    exp_beds = set([lib[x] for x in lib if re.match('exp_bed.*', x)])
    exp_beds = [x.rstrip('.bed') + '.txt' for x in exp_beds]
    experimentals = [x for x in cluster_file_list if (
        re.search('exp', x) or (os.path.basename(x) in exp_beds))]
    output_filename = lib['permutation_peaks_dir'] + '/combined_exp.txt'
    table = get_combined_clusters_from_file_list(experimentals, exons_as_rows, output_filename=output_filename)
    verification(table)
    control_beds = set([lib[x] for x in lib if re.match('control_bed.*', x)])
    control_beds = [x.rstrip('.bed') + '.txt' for x in control_beds]
    controls = [x for x in cluster_file_list if (
        re.search('control', x) or (os.path.basename(x) in control_beds))]
    output_filename = lib['permutation_peaks_dir'] + '/combined_control.txt'
    table = get_combined_clusters_from_file_list(controls, exons_as_rows, output_filename=output_filename)
    verification(table)
