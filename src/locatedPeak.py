import HTSeq
import re
import sys
import os
import pandas
import numpy as np


class locatedPeak(object):

    def __init__(self, chrm, left, right, strand, gene_name):
        self.left = left
        self.right = right
        self.strand = strand
        self.gene_name = gene_name
        self.chrm = chrm

    def iv_overlap(self, a, b):
        tup = self.overlap(a, b)
        _ivs = [None, None, None]
        for x in range(len(tup)):
            if tup[x] is not None:
                _ivs[x] = HTSeq.GenomicInterval(
                    a.chrom, tup[x][0], tup[x][1], a.strand)
        return _ivs

    def overlap(self, a, b):
        if a.start < b.start:
            if a.end < b.start:
#'''       b1--------b2
#   a1--a2'''
                return ((a.start, a.end), None, None)
            else:
                if a.end < b.end:
#'''       b1--------b2
#       a1-----a2'''
                    return ((a.start, b.start), (b.start, a.end), None)
                else:
#'''       b1--------b2
#       a1----------------a2'''
                    return ((a.start, b.start), (b.start, b.end),
                            (b.end, a.end))
        else:
            if a.start < b.end:
                if a.end < b.end:       
#'''       b1--------b2
#             a1--a2       '''
                    return (None, (a.start, a.end), None)
                else:
#'''       b1--------b2
#                a1-----a2'''
                    return (None, (a.start, b.end), (b.end, a.end))
            else:
#'''       b1--------b2
#                          a1-----a2'''
                return (None, None, (a.start, a.end))

                
    def add_location_from_integral(self, gtf, ga):
        left = self.left
        right = self.right
        strand = self.strand
        chrm = self.chrm
        print ">%s" % self.gene_name
        rows = gtf[self.gene_name]
        if len([1 for x in rows if x['2']=='CDS']) == 0:
            self.location = 'ncRNA'
            print rows
            return self.location
        cds_left = min([x['3'] for x in rows if x['2']=='CDS'])
        cds_right = max([x['4'] for x in rows if x['2']=='CDS'])
        cds_iv = HTSeq.GenomicInterval(
            self.chrm, cds_left, cds_right, self.strand)
        peak_iv = HTSeq.GenomicInterval(
            self.chrm, self.left, self.right, self.strand)
        overlaps = self.iv_overlap(peak_iv, cds_iv)
        integrals = {}
        for label, _iv in zip(['left_utr', 'cds', 'right_utr'], overlaps):
            if _iv is not None:
                integrals[label] = self.integral(ga, _iv)
            else:
                integrals[label] = 0
        print integrals
        location_l = sorted(zip(integrals.keys(), integrals.values()),
                          key=lambda x: x[1])
        print location_l
        location = location_l[-1][0]
        if strand == '+':
            if location == 'left_utr': self.location = "5'UTR"
            if location == 'right_utr': self.location = "3'UTR"
        if strand == '-':
            if location == 'left_utr': self.location = "3'UTR"
            if location == 'right_utr': self.location = "5'UTR"
        if location == 'cds': self.location = "CDS"
        print "*"
        print location, strand
        print self.location
        return self.location

    def integral(self, ga, iv):
        auc = 0
        try:
            for _iv, score in ga[iv].steps():
                auc += score * (_iv.end - _iv.start)
        except:  # Zero length interval. Not sure how to check this?
            return 0
        return auc
'''

--------------a==========b-------c========d------
----------x===z=======y--------------------------
for d in [a, b, c, d]:
    if x >= d: range = 0
    range = x -> min(d, y)
    integral = ga[range].integral
'''
