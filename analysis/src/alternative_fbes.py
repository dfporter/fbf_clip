import HTSeq
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
#from scipy import signal
import numpy as np
import random
import pickle
import pandas
import re
import os
import sys
import glob
#import pysam
#import scipy as sp
#import statsmodels.api as sm


def get_sequences(combined):
    fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for index, peak_row in combined.iterrows():
        start = combined.loc[index, 'left']
        end = combined.loc[index, 'right']
        chrm = combined.loc[index, 'chrm']
        seq = sequences[chrm][start:end]
        #print "%s:%i-%i: seq %s" % (chrm, start, end, seq)
        if combined.loc[index, 'strand'] == '-':
            seq = rc(seq)
        combined.loc[index, 'seq'] = seq


def get_peaks(filename):
    peaks = pandas.DataFrame.from_csv(filename, sep='\t')
    return peaks


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

def initialize_counts(range_len):
    upstream_counts = {}
    downstream_counts = {}
    for pos in range(0, 4):
        upstream_counts[pos] = {}
        for nt in ['n', 'a', 't', 'g', 'c']:
            upstream_counts[pos][nt] = 0
    for pos in range(0, range_len):
        downstream_counts[pos] = {}
        for nt in ['n', 'a', 't', 'g', 'c']:
            downstream_counts[pos][nt] = 0
    return (upstream_counts, downstream_counts)

def convert_counts_to_dict(upstream_seqs, downstream_seqs, range_len):
    upstream_counts, downstream_counts = initialize_counts(range_len)
    neg_upstream_counts, neg_downstream_counts = initialize_counts(range_len)
    for seq in upstream_seqs:
        for pos in range(0,4):
            upstream_counts[pos][seq[pos]] += 1
    for seq in downstream_seqs:
        for pos in range(0,range_len):
            downstream_counts[pos][seq[pos]] += 1
    return (upstream_counts, downstream_counts)

def get_pos_and_neg(peaks, out_dir, pos='', neg='', label=''):
    range_len = 3
    print "Target seq: %s. Control: %s." % (pos, neg)
    (upstream_seqs, downstream_seqs, fasta_lines) = get_targ_seq(peaks, pos, range_len)
    (neg_upstream_seqs, neg_downstream_seqs, neg_fasta_lines
     ) = get_targ_seq(peaks, neg, range_len)
    with open('%s/%s_%s.fa' % (out_dir, pos, label), 'w') as f:
        f.write(fasta_lines)
    with open('%s/%s_%s.fa' % (out_dir, neg, label), 'w') as f:
        f.write(neg_fasta_lines)
    upstream_counts, downstream_counts = convert_counts_to_dict(
        upstream_seqs, downstream_seqs, range_len)
    neg_upstream_counts, neg_downstream_counts = convert_counts_to_dict(
        neg_upstream_seqs, neg_downstream_seqs, range_len)
    print "upstream:"
    print upstream_counts
    print "downstream:"
    print downstream_counts
    upstream_fractions = convert_to_fractions(upstream_counts)
    downstream_fractions = convert_to_fractions(downstream_counts)
    neg_upstream_fractions = convert_to_fractions(neg_upstream_counts)
    neg_downstream_fractions = convert_to_fractions(neg_downstream_counts)
    print "upstream:"
    for pos in upstream_fractions:
        for nt in upstream_fractions[pos]:
            print "%i/%s: %s. vs neg: %f" % (
                pos, nt, str(upstream_fractions[pos][nt]),
                float(upstream_fractions[pos][nt]
                      )/float(neg_upstream_fractions[pos][nt]))
    print "downstream:"
    for pos in downstream_fractions:
        for nt in downstream_fractions[pos]:
            if neg_downstream_fractions[pos][nt] > 0:
                print "%i/%s: %s. vs neg: %f" % (
                    pos, nt, str(downstream_fractions[pos][nt]),
                    float(downstream_fractions[pos][nt]
                          )/float(neg_downstream_fractions[pos][nt]))
            else:
                print "%i/%s: %s. no neg value (%s)" % (
                    pos, nt, str(downstream_fractions[pos][nt]),
                    str(neg_downstream_fractions[pos][nt]))

def get_targ_seq(peaks, target_seq, range_len):
    fasta_lines = ""
    num_motif = 0
    upstream_seqs = []
    downstream_seqs = []
    intermediate_len = range_len
    for index, row in peaks.iterrows():
        seq = peaks.loc[index, 'seq']
        m = re.finditer(
            '(.{0,4})' + target_seq + '(.{1,%i})at' % int(intermediate_len),
                        seq)
        for instance in m:
            if len(instance.group(1)) == 4:
                upstream = instance.group(1)
            else:
                upstream = 'n' * (4 - len(instance.group(1)))
                upstream += instance.group(1)
            upstream_seqs.append(upstream)
            if instance.group(2) == intermediate_len:
                downstream = instance.group(2)
            else:
                downstream = 'n' * (intermediate_len - len(instance.group(2)))
                downstream += instance.group(2)
            downstream_seqs.append(downstream)
            num_motif += 1
            fasta_lines += ">%s|%s\n%s\n" % (
                str(num_motif),
                peaks.loc[index, 'gene_name'],
                instance.group(0))
    return (upstream_seqs, downstream_seqs, fasta_lines)

def convert_to_fractions(counts):
    totals = {}
    fractions = {}
    for pos in counts:
        fractions[pos] = {}
        totals[pos] = sum([counts[pos][x] for x in counts[pos]])
        for nt in counts[pos]:
            fractions[pos][nt] = float(counts[pos][nt])/float(totals[pos])
    return fractions

def overlap_fbe(fbe_pos, m):
    for a_fbe_pos in fbe_pos:
        if a_fbe_pos[0] <= m.start() <= a_fbe_pos[1]:
            return True
        if a_fbe_pos[0] <= m.end() <= a_fbe_pos[1]:
            return True
    return False


def other_sites(peaks, out_dir, num_shuffles=10, exclude_fbes=True):
    counts = {}
    rand_counts = {}
    seqrs_motifs = ['gtgt\w\w\wtg', 'gtgt\w\w\wtt', 'tgtgtatata',
              'gttgt\w\w\wtg',
              'gttgt\w\w\wtt',
              'tgt\w\w\wat',]
    seqrs_motifs_keys = ['gtgtnnntg', 'gtgtnnntt', 'tgtgtatata',
              'gttgtnnntg',
              'gttgtnnntt',
              'tgtnnnat',]
    motifs = ['tgt', 'ctgt', 'c\wtgt',
              '[^c]ctgt', 'c[^c]tgt',
              'c\w?tgt',
              'ttg', 'gtt', 'cttg', 'cgtt']
    motifs.extend(seqrs_motifs)
    fbe = 'tgt[^c][^g][^g]at'
    motifs_keys = ['tgt', 'ctgt', 'cntgt',
                  'only -1C TGT',
                  'only -2C TGT',
                  '-1 or -2C TGT',
                  'ttg', 'gtt', 'cttg', 'cgtt']
    motifs_keys.extend(seqrs_motifs_keys)
    # Oops.
    for motif in motifs_keys:
        counts[motif] = {'with': 0, 'without': 0}
        for index in range(0, num_shuffles):
            rand_counts.setdefault(index, {})
            rand_counts[index][motif] = {'with': 0, 'without': 0}
    for index, row in peaks.iterrows():
        if index % 100 == 0:
            print index
        seq = peaks.loc[index, 'seq'].lower()
        fbe_pos = []
        for m_fbe in re.finditer(fbe, seq):
            fbe_pos.append((m_fbe.start(), m_fbe.end()))
        for motif in motifs_keys:
            motif_regex = motifs[motifs_keys.index(motif)]
            has_match = False
            for m in re.finditer(motif_regex, seq):
                if exclude_fbes and overlap_fbe(fbe_pos, m):
                    continue
                has_match = True
            if has_match:
                counts[motif]['with'] += 1
            else:
                counts[motif]['without'] += 1
            for index in range(0, num_shuffles):
                seqarr = list(seq)
                random.shuffle(seqarr)
                rand_seq = ''.join(seqarr)
                has_match = False
                for m2 in re.finditer(motif_regex, rand_seq):
                    if exclude_fbes and overlap_fbe(fbe_pos, m):
                        continue
                    has_match = True
                if has_match:
                    rand_counts[index][motif]['with'] += 1
                else:
                    rand_counts[index][motif]['without'] += 1
    for motif in motifs_keys:
        num_above = 0
        num_below = 0
        print motif
        print 'exp %s' % str(counts[motif])
        for index in rand_counts:
            print 'rand %i' % rand_counts[index][motif]['with']
            if rand_counts[index][motif]['with'] >= counts[motif]['with']:
                num_above += 1
#                print 'increment above'
            else:
                num_below += 1
#                print 'increment below'
        counts[motif]['num_above'] = num_above
        counts[motif]['num_below'] = num_below
        counts[motif]['pvalue'] = float(num_above)/float(
            num_above + num_below)
        print counts[motif]
    return (counts, rand_counts)


def write_fastas(peaks, output_basename):
    li = ''
    no_fbe_li = ''
    for index, row in peaks.iterrows():
        fasta_line = '>{n}\n{s}\n'.format(
            n=str(str(index) + '_' + peaks.loc[index, 'gene_name']),
            s=peaks.loc[index, 'seq'])
        li += fasta_line
        if re.search('tgt\w\w\wat', peaks.loc[index, 'seq'].lower()) is None:
            no_fbe_li += fasta_line
    if not os.path.exists('data/'): os.system('mkdir data/')
    if not os.path.exists('data/fasta'): os.system('mkdir data/fasta')
    with open('data/fasta/' + output_basename, 'w') as f: f.write(li)
    if not os.path.exists('data/fasta/no_fbes/'):
        os.system('mkdir data/fasta/no_fbes')
    with open('data/fasta/no_fbes/' + output_basename, 'w') as f:
        f.write(no_fbe_li)

from add_biotype import *

if __name__ == '__main__':
    num_shuffles = 100
    top_level_dir = sys.argv[1] + '/'
    out_dir = 'data/' + os.path.basename(top_level_dir.rstrip('/')) + '/'
    peaks = {}
    for filename in glob.glob(top_level_dir + 'combined*.txt'):
        peaks[filename] = pandas.read_csv(filename, sep='\t')
        if 'biotype' not in peaks[filename].columns:
            add_biotype_to_all(glob.glob(top_level_dir + 'combined*.txt'))
            break
    for filename in glob.glob(top_level_dir + 'combined*.txt'):
        print "Loading %s..." % filename
        peaks[filename] = pandas.read_csv(filename, sep='\t')
        if 'biotype' not in peaks[filename].columns:
            from add_biotype import *
            add_biotype_to_all(glob.glob(top_level_dir + 'combined*.txt'))
            peaks[filename] = pandas.read_csv(filename, sep='\t')
        peaks[filename] = peaks[filename][peaks[filename]['biotype']=='protein_coding']
        get_sequences(peaks[filename])
        out_basename = os.path.basename(filename).partition('.txt')[0] + '.fa'
        write_fastas(peaks[filename], out_basename)
    #sys.exit()
    for filename in peaks:
        print "Examining %s..." % filename
        if not os.path.exists('tables/{d}'.format(
            d=os.path.basename(os.path.dirname(filename))
        )):
            os.system('mkdir tables/{d}'.format(
                d=os.path.basename(os.path.dirname(filename))
            ))
        out_f = open('tables/{d}/cryptic_fbes_{num_shuf}_{b}'.format(
            d=os.path.basename(os.path.dirname(filename)),
            b=os.path.basename(filename),
            num_shuf=num_shuffles), 'w')

        #get_pos_and_neg(peaks[filename], out_dir, pos='tgt',
        #         neg='ttg', label=os.path.basename(filename))
        (counts, rand_counts) = other_sites(peaks[filename], out_dir,
                                            num_shuffles=num_shuffles)
        li = "%s:\n" % filename
        for motif in counts:
            sum_motif = float(counts[motif]['with'] + counts[motif]['without'])
            li += "\t%s:\n" % motif
            li += "\t\twith %f (%i) without %f (%i)\n" % (
                float(counts[motif]['with'])/float(sum_motif),
                counts[motif]['with'],
                float(counts[motif]['without'])/float(sum_motif),
                counts[motif]['without'])
            li += "\t\trandom peak seq: with %f (%i) without %f (%i)\n" % (
                float(rand_counts[0][motif]['with'])/float(sum_motif),
                rand_counts[0][motif]['with'],
                float(rand_counts[0][motif]['without'])/float(sum_motif),
                rand_counts[0][motif]['without'])
            if rand_counts[0][motif]['with'] > 0:
                li += "\t\tenrichment: %f. p value from %i simulations: %f\n" % (
                    float(counts[motif]['with'])/float(
                        rand_counts[0][motif]['with']),
                    num_shuffles,
                    float(counts[motif]['pvalue']))
            else: li += "\t\tNo counts in negative.\n"
        out_f.write(li)
        out_f.close()
