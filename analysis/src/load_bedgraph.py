import HTSeq
import os
import collections
import re


def load_bedgraph(fname):
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    plus_file = fname.partition('.wig')[0] + '_+.wig'
    add_strand_to_ga_from_bedgraph_file(plus_file, ga, '+')
    minus_file = fname.partition('.wig')[0] + '_-.wig'
    add_strand_to_ga_from_bedgraph_file(minus_file, ga, '-')
    return ga


def add_strand_to_ga_from_bedgraph_file(fname, ga, strand):
    with open(fname, 'r') as f:
        next(f)
        for line in f:
            #if line[0] != 'M': continue
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), strand)] = float(s[3])
    return ga


def get_val(ga, iv):
    return np.max(np.fromiter(ga[iv], dtype=np.float))


def get_bed_size(fname):
    return float(len(open(fname).readlines()))

