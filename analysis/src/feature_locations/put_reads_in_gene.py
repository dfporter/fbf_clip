import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys
import numpy as np


def assign_reads_to_genes(ga, gtf_fname):
    gtf = HTSeq.GFF_Reader(gtf_fname, end_included=True)
    gene_ids = set()
    reads_in_gene = {}
    for fname in ga:
        print fname
        for iv, val in ga[fname].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            reads_in_gene.setdefault(gene_id, ReadsInGene_obj(gene_id))
            reads_in_gene[-1].reads.append((iv.start, iv.end))
    stats_of_readsInGene_list(reads_in_gene.values())
    return reads_in_gene


def stats_of_readsInGene_list(_list):
    num_reads = {}
    for rig in _list:
        num_reads[rig.gene_id] = len(rig.reads)
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


class ReadsInGene:

    def __init__(self, gene_id):
        self.gene_id = gene_id
        self.reads = []

    def define_extent(self):
        pass

    def convert_reads_to_relative_location(self):
        pass

    def shuffle_reads_to_call_peaks(self):
        pass
