import HTSeq
import re
import sys
import os
import pandas
import numpy as np

from peaks import peaks
from locatedPeak import locatedPeak
sys.path.insert(0,
                '/groups/Kimble/Aman Prasad/clip2/analysis/src/')
import add_info_columns

def seq_from_iv(chrm, start, end, strand, sequences):
    seq = sequences[chrm][start:end]
    if strand == '-':
        return rc(seq)
    else:
        return seq

    
def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s


class annotatedPeaks(peaks):

    def add_biotype(self, gtf):
        to_biotype = {} #collections.defaultdict(set)
        for tup in zip(gtf.gene_name, gtf.biotype):
            if (tup[0] in to_biotype) and (
                to_biotype[tup[0]] == 'protein_coding'):
                continue
            to_biotype[tup[0]] = tup[1]
        self.data['biotype'] = [to_biotype[x] for x in self.data.gene_name]
        
    def add_seqs(self, sequences):
        seq = [seq_from_iv(tup[0], tup[1], tup[2], tup[3], sequences) \
               for tup in \
         zip(self.data.chrm, self.data.left, self.data.right, self.data.strand)]
        self.data['seq'] = seq

    def add_fbe(self, sequences=None):
        if (not hasattr(self, 'seq')) and (sequences is not None):
            self.add_seqs(sequences)
        pat = re.compile('tgt\w\w\wat', re.IGNORECASE)
        self.data['has_fbe'] = [has_pat(x, pat) for x in self.data.seq]
        self.data['number_of_fbes_fbe'] = [num_pat(x, pat) for x in self.data.seq]
        pat = re.compile('ctgt\w\w\wat', re.IGNORECASE)
        self.data['minus_one_c'] = [has_pat(x, pat) for x in self.data.seq]
        pat = re.compile('c\wtgt\w\w\wat', re.IGNORECASE)
        self.data['minus_two_c'] = [has_pat(x, pat) for x in self.data.seq]
        self.data['minus_one_or_two_c'] = [either(tup[0], tup[1]) \
            for tup in zip(self.data.minus_one_c, self.data.minus_two_c)]

    def add_location(self, gtf, use_this_column=None):
        add_info_columns.locate_in_gene(gtf, self.data)

    def add_location_from_integral(self, gtf, ga):
        ivs = zip(self.data.chrm, self.data.left, self.data.right,
                  self.data.strand, self.data.gene_name)
        located_peaks = [locatedPeak(*iv) for iv in ivs]
        locations = [x.add_location_from_integral(gtf, ga) for \
                     x in located_peaks]
        self.data['location'] = locations
        
    def copy_columns(self):
        self.data['n2_from_fbf1_reads'] = self.data['fbf1_n2']
        self.data['n2_from_fbf2_reads'] = self.data['fbf2_n2']
        self.data['fbf1_n2_reads'] = self.data['fbf1_n2']
        self.data['fbf2_n2_reads'] = self.data['fbf2_n2']
        self.data['rna_seq_reads'] = self.data['rna_seq']

    def write(self, filename):
        self.data.to_csv(filename, sep='\t', index=False)

def has_pat(seq, pat):
    if pat.search(seq.lower()):
        return 1
    else:
        return 0

def num_pat(seq, pat):
    if pat.search(seq.lower()):
        return len(re.findall(pat, seq.lower()))
    else:
        return 0

def either(x, y):
    if str(x) == '0' and str(y) == '0': return 0
    return 1
