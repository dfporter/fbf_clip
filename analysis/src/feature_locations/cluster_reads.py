import collections
import sys
import HTSeq
import numpy as np
import bisect
import pandas
import copy
import os
from statsmodels.sandbox.stats.multicomp import multipletests
import p_values_of_clusters


def is_good_gene(rows):
    # print rows
    for row in rows:
        if row['biotype'] == 'rRNA':
            return False
        if row['biotype'] == 'snoRNA':
            return False
        # if row['biotype'] == 'protein_coding':
        #     return False
        # if row['biotype'] != 'protein_coding':
        # # if row['gene_id'] != 'WBGene00010957':
        #     return False
    return True


def run(reads_by_gene_by_fname, counts_by_gene, exons, lib=None, args=None):
    """
    Find overlapping reads in a given gene from a list of ivs.
Input:
    reads_by_gene_by_fname: A dict by filename, then by gene_id, holding a list of ivs.
Internal:
    rigs: A list of iv objects, one for each read. Includes all reads in a given
            gene.
    clusters: A list of [start, stop, max read overlap] for each region of
            continuous coverage.
    txpt_exons: [(exon start, exon stop)...] for one txpt.
    clusters_w_pvalue: A dict with key=gene_id and value= A list of
            [start, stop, max read overlap, p value] for each region of
            continuous coverage.
Output:
    table: Dataframe of peaks.
    """
    exons = dict((k, v) for k, v in exons.iteritems() if is_good_gene(v))
    for fname in reads_by_gene_by_fname:
        print "*" * 7 + '\n' + fname
        reads_by_gene = reads_by_gene_by_fname[fname]
        stats(reads_by_gene)
        print str(reads_by_gene.keys())[:100]
        print str(exons.keys())[:100]
        gene_list = set(set(reads_by_gene.keys()) & set(exons.keys()))
        # gene_list = ['WBGene00001595', 'WBGene00014568']
        if len(gene_list) == 0:
            print "No targets found for {i}. Skipping...".format(i=fname)
            continue
        get_p_values_from_scramble_function = p_values_of_clusters.get_exon_ranges_and_scramble_and_return_probabilities
        clusters_w_pvalue = dict(zip(
            gene_list,
            [add_p_values(*tup) for tup in yield_gene_with_tup(
                gene_list, reads_by_gene, exons, get_p_values_from_scramble_function)]
        ))
        # for gene in clusters_w_pvalue:
        #     print "{g}: {v}".format(g=gene, v=clusters_w_pvalue[gene])
        table = convert_to_table(clusters_w_pvalue, exons)
        # Get rid of clusters with only one read.
        table = table[table['max_coverage'] > 1]
        if len(table.index) == 0:
            print 'empty table'
            continue
        add_stats(table)
        output_filename = lib['clusters_dir'] + '/' + os.path.basename(
            fname).partition('.bed')[0] + '.txt'
        write_table(table, output_filename=output_filename)
        print table.head(3)
    return table


def write_table(table, output_filename='temp.txt'):
    if not os.path.exists(os.path.dirname(output_filename)):
        os.system('mkdir ' + os.path.dirname(output_filename))
    expected_columns = [
        'chrom', 'left', 'right', 'max_coverage', 'gene_id', 'strand', 'pval',
    ]

    extra_columns = list(set(table.columns) - set(expected_columns))
    column_order = expected_columns + extra_columns
    table = table[column_order]
    table.to_csv(output_filename, sep='\t', index=False)


def convert_to_table(dict_by_gene_of_clusters, exons):
    expected_columns = [
        'chrom', 'left', 'right', 'max_coverage', 'gene_id', 'strand', 'pval'
    ]
    table = []
    for gene in dict_by_gene_of_clusters:
        if isinstance(dict_by_gene_of_clusters[gene][0], type('a')):
            # Ignore NA.
            continue
        for row in dict_by_gene_of_clusters[gene]:
            if row is None: continue
            if row[0] != 'NA':
                table.append([gene, row[0], row[1], row[2], row[3]])
    if len(table) == 0:
        df = pandas.DataFrame(columns=expected_columns)
        return df
    table = [
        {'chrom': exons[gene][0]['0'], 'strand': exons[gene][0]['6'],
         'gene_id': gene, 'left': left, 'right': right, 'max_coverage': coverage,
         'pval': pval} for (gene, left, right, coverage, pval) in table
    ]
    table = pandas.DataFrame(table)
    try:
        table = table[expected_columns]
    except:
        print "cluster_reads.convert_to_table(): Error?"
        print table
    return table


def add_stats(table, alpha=0.01):
    past_fdr, padj, _, _ = multipletests(table['pval'].astype(float), alpha=alpha, method='fdr_bh')
    table['padj'] = padj
    table['past_fdr'] = past_fdr


def yield_gene_with_tup(genes, rbg, exons, get_p_values_from_scramble_function):
    for gene in genes:
        yield gene, rbg, exons, get_p_values_from_scramble_function


def add_p_values(gene, reads_by_gene, exons, get_pvals_from_scramble):
    rig = reads_by_gene[gene]
    if len(rig) < 2:
        return ['NA', 'NA', '0.0', '1.0']
    clusters = p_values_of_clusters.find_clusters_in_gene(rig)
    if len(clusters) < 1:
        return ['NA', 'NA', '0.0', '1.0']
    if max([x[2] for x in clusters]) < 2:
        return ['NA', 'NA', '0.0', '1.0']
    txpt_exons = exons[gene]
    txpt_id, txpt_exons = p_values_of_clusters.find_longest_txpt(txpt_exons)
    probabilities = get_pvals_from_scramble(
        num_reads=len(rig), num_permutations=1000,
        txpt_id=txpt_id, exons=txpt_exons)
    clusters_w_pvalue = p_values_of_clusters.add_p_values_from_histogram(
        probabilities, clusters
    )
    return clusters_w_pvalue


def stats(rbg):
    if len(rbg.keys()) == 0: return
    num_reads = {}
    for gene in rbg:
        num_reads[gene] = len(rbg[gene])
    total_reads = sum([num_reads[x] for x in num_reads])
    total_genes = len(list(num_reads.keys()))
    print """
    Genes: {n}
    Reads: {r}
    Mean reads per gene: {mr}
    Median reads per gene: {medr}
    Highest number per gene: {h}
    Lowest number per gene: {s}
    """.format(
        n=total_genes, r=total_reads,
        mr=np.mean(num_reads.values()),
        medr=np.median(num_reads.values()),
        h=max(num_reads.values()),
        s=min(num_reads.values()),
    )

