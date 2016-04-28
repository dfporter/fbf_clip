__author__ = 'dp'
import re
import os
import sys


fastx_collapser_fasta_file = sys.argv[1]
fastq_file = sys.argv[2]


seq_dups = set()
with open(fastx_collapser_fasta_file, 'r') as f:
    while True:
        id = f.readline()
        seq = f.readline()
        if not seq: break
        seq_dups.add(seq)
outf = 'fitered.fastq'
added_dups = set()
with open(fastq_file, 'r') as f:
    read_name = f.readline()
    seq = f.readline()
    _ = f.readline()
    qual = f.readline()
    if seq not in added_dups:
        if seq in seq_dups:
            outf.write("".join([read_name, seq, _, qual]))
            seq_dups.remove(seq)
            added_dups.add(seq)
        else:
            outf.write("".join([read_name, seq, _, qual]))

