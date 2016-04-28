import HTSeq
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import pickle
import pandas
import re
import os
import glob
import pysam


def get_sequences():
    fasta_filename = '/home/dp/Desktop/ensemblEF4/ef4.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    return sequences

def create_gtf_with_names_file(gtf_noheader_filename):
    outf = open('gtf_with_names_column.txt', 'w')
    with open(gtf_noheader_filename, 'r') as f:
        outf.write('\t'.join([str(x) for x in range(0,9)]) + '\tgene_name\ttranscript_id\ttranscript_name\tgene_id\texon_number\n')
        for li in f:
            s = li.rstrip('\n').split('\t')
            s = s + [get_name(li), get_txpt_name(li), get_txpt_id(li), get_gene_id(li), get_exon_num(li)]
            outf.write("\t".join(s) + '\n')
    outf.close()
    gtf_df = pandas.read_csv('gtf_with_names_column.txt', sep='\t', header=0)
    return gtf_df
def get_name(line):
    #line = str(line)
    m = re.search('gene_name "([^"]+)"', line)
    if m is not None:
        return m.group(1)
    else:
        return 'NA'
def get_txpt_id(line):
    #line = str(line)
    m = re.search('transcript_id "([^"]+)"', line)
    if m is not None:
        return m.group(1)
    else:
        return 'NA'
def get_gene_id(line):
    #line = str(line)
    m = re.search('gene_id "([^"]+)"', line)
    if m is not None:
        return m.group(1)
    else:
        return 'NA'
def get_txpt_name(line):
    #line = str(line)
    m = re.search('transcript_name "([^"]+)"', line)
    if m is not None:
        return m.group(1)
    else:
        return 'NA'
def get_exon_num(line):
    #line = str(line)
    m = re.search('exon_number "([^"]+)"', line)
    if m is not None:
        return m.group(1)
    else:
        return 'NA'
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
def seq_from_iv(iv, sequences):
    """Returns a slice from 0-based, slice coordinates.
    """
    if iv[2] > len(sequences[iv[0]]):
        end = len(sequences[iv[0]])
    else:
        end = iv[2]
    a_seq = sequences[iv[0]][iv[1]:end]
    if iv[3] == "-":
        a_seq = rc(a_seq)
    return a_seq

#def gene:
#    def __init__(self, exons):
#        self.exons_list = exons
def has_intron(txpt_id, gtf_df):
    if len(gtf_df[(gtf_df['transcript_id']==txpt_id) & (gtf_df['2']=='CDS')]) > 1:
        return True
    return False


#gtf_noheader_filename = "/home/dp/Desktop/celegans_genome/wormbase_ws235/Caenorhabditis_elegans.WBcel235.78.noheader.gtf" 
gtf_noheader_filename = "/home/dp/Desktop/ensemblEF4/Saccharomyces_cerevisiae.EF4.70.gtf"
gtf_df = create_gtf_with_names_file(gtf_noheader_filename)
seqs = get_sequences()

w_df = pandas.DataFrame.from_csv('../clip_results/src/predictTargets/predictTargets/lib/weissman.txt', sep='\t')
names = set(w_df['common.name'])
gtf_df['in_weis'] = 0
gtf_df_full = gtf_df.copy()
for index, row in gtf_df.iterrows():
    if str(gtf_df.loc[index, 'gene_name']) in names:
        gtf_df.loc[index, 'in_weis'] = 1
gtf_df = gtf_df[gtf_df['in_weis']==1]
txpt_ids = set(gtf_df['transcript_id'])
txpt_ids = [x for x in txpt_ids if not pandas.isnull(x)]
seqs = get_sequences()
seqs.keys()
txpt_seqs, upstream, downstream, combined = add_gene_seq(gtf_df, seqs)
import regex

num_motifs = {'TAAT': [], 'TGTANATA': []}
for gene in txpt_seqs:
    num_motifs_puf2 = len(regex.findall('TAAT', combined[gene], overlapped=True))
    num_motifs_classical = len(regex.findall('TGTA\wATA', combined[gene], overlapped=True))
    num_motifs['TAAT'].append(num_motifs_puf2)
    num_motifs['TGTANATA'].append(num_motifs_classical)

import scipy as sp
import numpy as np
arr = np.asarray(num_motifs['TAAT'])
print sp.mean(arr)
print sp.median(arr)

lengths = []
for gene in txpt_seqs:
    lengths.append(len(combined[gene]))
lengths = np.asarray(lengths)
print sp.mean(lengths)
print sp.median(lengths)

def signal_at_pos(gene, gtf_df, rel_pos, coverage, sequences, use_bam):
    """Bam files and our sequence string are 0-based.
    gtf_df (.gtf format) is 1-based.
    """
    stop_codon_df = gtf_df[(gtf_df['gene_name']==gene) & (gtf_df['2']=='stop_codon')]
    indx = stop_codon_df.index[0]
    if len(stop_codon_df) == 0:
        print "No stop codon for %s." % gene
        return 'NA'
    else:
        strand = gtf_df.loc[indx, '6']
        if strand == '+':
            # chr_stop = int(gtf_df.loc[indx, '3']) - 1
            # slice of stop codon is sequences[chr_stop:chr_stop+3]
            chr_pos = int(gtf_df.loc[indx, '3']) - 2 + rel_pos
            # slice of genomic sequence is chr_pos:+4 is TAAT.
            # Not sure why there is an offset here, but it slices correctly.
            iv = (gtf_df.loc[indx, '0'], chr_pos, chr_pos + 4, strand)
            #if sequences[gtf_df.loc[indx, '0']][chr_pos:chr_pos+4] != 'TAAT':
            #    logfile = open('logfile', 'a')
            #    error_line = "Error %s +, rel_pos %i not a TAAT %s." % (gene, rel_pos, sequences[gtf_df.loc[indx, '0']][chr_pos-1:chr_pos+5])
            #    print error_line
            #    logfile.write(error_line + '\n')
            #    logfile.close()
        else:
            chr_pos = int(gtf_df.loc[indx, '3']) - 1 - rel_pos  # -1 to change from 1-based .gtf.
            # slice of genomic sequence = chr_pos:chr_pos+4 is TAAT.
            #if sequences[gtf_df.loc[indx, '0']][chr_pos:chr_pos+4] != 'ATTA':
            #    logfile = open('logfile', 'a')
            #    error_line = "Error %s -, rel_pos %i not a TAAT %s." % (gene, rel_pos, sequences[gtf_df.loc[indx, '0']][chr_pos-1:chr_pos+5])
            #    print error_line
            #    logfile.write(error_line + '\n')
            #    logfile.close()
            iv = (gtf_df.loc[indx, '0'], chr_pos, chr_pos+4, strand)
        return signal_at_iv(iv, coverage, use_bam)

def signal_at_iv(iv, coverage, use_bam):
    """If use_bam is False, coverage is an HTSeq genomic array.
    If use_bam is True, coverage is a bam filename.
    """
    if not use_bam:
        window = HTSeq.GenomicInterval(iv[0], max(0, iv[1]), iv[2], iv[3])
        return max(np.fromiter(coverage[window], dtype='i'))
    left_boundary = iv[1]
    right_boundary = iv[2]
    chrom = iv[0]
    strand = iv[3]
    bamfile = pysam.AlignmentFile(coverage, "rb")
    s = bamfile.fetch(chrom, left_boundary, right_boundary)
    total_reads = 0
    for r in s:
        if(strand=="+" and not r.is_reverse):
             total_reads += 1
        if(strand=="-" and r.is_reverse):
             total_reads += 1
    bamfile.close()
    return total_reads

def get_pos(gtf_df, re_motif, coverage, sequences, use_bam=False):
    """
    sequences are only used for debugging puproses (the final argument, I mean).
    """
    pos = []
    #gtf_df = gtf_df[gtf_df['0']=='I']
    value_at_pos = {}
    #for gene in combined:
    for gene in set(gtf_df['gene_name']):
        print gene
        if gene not in combined:
            print "Has introns, skipping..."
            continue
        if list(gtf_df[(gtf_df['gene_name']==gene)]['0'])[0]=='Mito':
            print "Skipping mtDNA gene..."
            continue
        a_seq = combined[gene]
        all_m = regex.finditer(re_motif, a_seq, overlapped=True)
        for motif in all_m:
            if strand == '-':
                stop_pos = len(a_seq) - 200
                rel_pos = motif.end() - stop_pos
                #substr = a_seq[rel_pos+stop_pos-4:rel_pos+stop_pos]  # This is TAAT.
                #print substr
                #stop_codon_seq = a_seq[stop_pos-3:stop_pos]
                pos.append(rel_pos)
                value_at_pos.setdefault(rel_pos, [])
                value_at_pos[rel_pos].append(signal_at_pos(gene, gtf_df, rel_pos, coverage, sequences, use_bam=use_bam))
            else:
                stop_pos = len(a_seq) - 200
                #stop_codon_seq = a_seq[stop_pos-3:stop_pos]
                rel_pos = motif.start() - stop_pos
                pos.append(rel_pos)
                value_at_pos.setdefault(rel_pos, [])
                value_at_pos[rel_pos].append(signal_at_pos(gene, gtf_df, rel_pos, coverage, sequences, use_bam=use_bam))
    return pos, value_at_pos
#gtf_all = gtf_df.copy()
def read_coverage(bam_filename):
    _bamfile = HTSeq.BAM_Reader(bam_filename)
    coverage = HTSeq.GenomicArray("auto", stranded=True, typecode='i')
    #gtf_df = pandas.read_csv(gtf_noheader_filename, sep='\t', header=None)
    print "Reading alignments from bamfile..."
    for aln in _bamfile:  # Very slow.
        if aln.aligned:
            coverage[aln.iv] += 1
    return coverage

def bamfile_size(bamfile_name):
        # One-liner from https://www.biostars.org/p/1890/:
    return reduce(lambda x, y: x + y,
            [eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamfile_name)])

bamfile = '/home/dp/Desktop/bams/berger_srr1066657.bam'
bamfile_sizes = {}
bamfile_sizes['rna_seq'] = bamfile_size(bamfile)
#rna_seq_cov = read_coverage(bamfile)
pos_taat, value_at_pos_taat = get_pos(gtf_df, 'TAAT', bamfile, sequences, use_bam=True)
pos_taag, value_at_pos_taag = get_pos(gtf_df, 'TAAG', bamfile, sequences, use_bam=True)
wt_bamfile = '/home/dp/Desktop/bams/polyN.bam'
wt_cov = read_coverage(wt_bamfile)
pos_taat_wt, value_at_pos_taat_wt = get_pos(gtf_df, 'TAAT', wt_cov, sequences)
r1_bamfile = '/home/dp/Desktop/bams/R1_SNE.bam'
r1_cov = read_coverage(r1_bamfile)
pos_taag_r1, value_at_pos_taag_r1 = get_pos(gtf_df, 'TAAG', r1_cov, sequences)
bamfile_sizes['polyN'] = bamfile_size(wt_bamfile)
bamfile_sizes['r1'] = bamfile_size(r1_bamfile)
normfactors = {
    'rna_seq': float(1e6 * 1.0/float(bamfile_sizes['rna_seq'])),
    'polyN': float(1e6 * 1.0/float(bamfile_sizes['polyN'])),
    'r1': float(1e6 * 1.0/float(bamfile_sizes['r1']))
    }

print "taat"
_x_taat, _y_taat = as_array(
    value_at_pos_taat, norm_factor=1, read_num_cutoff=2
    )
print "taag"
_x_taag, _y_taag = as_array(
    value_at_pos_taag, norm_factor=1, read_num_cutoff=2
    )
print "taat wt"
_x_taat_wt, _y_taat_wt = as_array(
    value_at_pos_taat_wt, norm_factor=1,
    )
print "taag r1"
_x_taag_r1, _y_taag_r1 = as_array(
    value_at_pos_taag_r1, norm_factor=1,
    )


def as_array(_value_at_pos, norm_factor=1.0, read_num_cutoff=-100):
    _x = []; _y = []; _x_near = []; _y_near = []
    n_excluded = 0
    n_input = 0
    for rel_pos in _value_at_pos:
        for read_num in _value_at_pos[rel_pos]:
            n_input += 1
            _x.append(rel_pos)
            if abs(rel_pos) < 500:
                if read_num > read_num_cutoff:
                    _x_near.append(rel_pos)
                    _y_near.append(norm_factor * read_num)
                else:
                    pass
            else:
                n_excluded += 1
            _y.append(read_num)
    print "n_input %i num_excluded by distance %i" % (n_input, n_excluded)
    x_arr = np.asarray(_x_near)
    y_arr = np.asarray(_y_near)
    return _x_near, _y_near


import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
#_x = np.linspace(-200, 200, num=401)
matplotlib.rcParams.update({'font.size': 8})
ax = plt.subplot(2,2,1)
plt.xlim([-500,200]); plt.ylim([-1.2,3])
plt.xlabel("Distance to a stop codon"); plt.ylabel("Site abundance")
plt.xticks((-400., -200., 0., 200.)); plt.yticks((-1.,0.,1.,2.))
plt.scatter(_x_taat, np.log10([_y * normfactors['rna_seq'] for _y in _y_taat]), alpha=0.01, edgecolor='none', color='black')
ax = plt.subplot(2,2,2)
plt.xlim([-500,200]); plt.ylim([-1.2,3])
plt.xlabel("Distance to a stop codon"); plt.ylabel("Site abundance")
plt.xticks((-400., -200., 0., 200.)); plt.yticks((-1.,0.,1.,2.,))
plt.scatter(_x_taag, np.log10([_y * normfactors['rna_seq'] for _y in _y_taag]), alpha=0.01, edgecolor='none', color='black')
ax = plt.subplot(2,2,3)
plt.xlim([-500,200]); plt.ylim([0,3])
plt.xlabel("Distance to a stop codon"); plt.ylabel("WT peak height")
plt.xticks((-400., -200., 0., 200.)); plt.yticks((0,1,2,3))
plt.scatter(_x_taat_wt, np.log10([_y * normfactors['polyN'] for _y in _y_taat_wt]), alpha=0.01, edgecolor='none', color='black')
ax = plt.subplot(2,2,4)
plt.xlim([-500,200]); plt.ylim([0,3])
plt.xlabel("Distance to a stop codon"); plt.ylabel("R1 SNE peak height")
plt.xticks((-400., -200., 0., 200.)); plt.yticks((0,1,2,3))
plt.scatter(_x_taag_r1, np.log10([_y * normfactors['r1'] for _y in _y_taag_r1]), alpha=0.01, edgecolor='none', color='black')
plt.show()

# Do with targets only.
targs_wt = pandas.DataFrame.from_csv(
    '/home/dp/Desktop/clip_results/exps/cerevisiae/polyN/peaksPastFdr.info',
    sep='\t')
targs_r1 = pandas.DataFrame.from_csv(
    '/home/dp/Desktop/clip_results/exps/cerevisiae/R1_SNE/peaksPastFdr.info',
    sep='\t')
targs_wt_list = set(targs_wt['gene'])
targs_r1_list = set(targs_r1['gene'])
gtf_df['wt_targ'] = 0
gtf_df['r1_targ'] = 0
gtf_df['wt_targ'] = gtf_df.apply(lambda row: row['gene_id'] in targs_wt_list, axis=1)
gtf_df['r1_targ'] = gtf_df.apply(lambda row: row['gene_id'] in targs_r1_list, axis=1)
gtf_wt_df = gtf_df[gtf_df['wt_targ']==1]
gtf_r1_df = gtf_df[gtf_df['r1_targ']==1]

bamfile = '/home/dp/Desktop/bams/berger_srr1066657.bam'
pos_taat_targs, value_at_pos_taat_targs = get_pos(gtf_wt_df, 'TAAT', bamfile)
pos_taag_targs, value_at_pos_taag_targs = get_pos(gtf_r1_df, 'TAAG', bamfile)
wt_bamfile = '/home/dp/Desktop/bams/polyN.bam'
pos_taat_wt_targs, value_at_pos_taat_wt_targs = get_pos(gtf_wt_df, 'TAAT', wt_bamfile)
r1_bamfile = '/home/dp/Desktop/bams/R1_SNE.bam'
pos_taag_r1_targs, value_at_pos_taag_r1_targs = get_pos(gtf_r1_df, 'TAAG', r1_bamfile)

_x_taat_targs, _y_taat_targs = as_array(value_at_pos_taat_targs)
_x_taag_targs, _y_taag_targs = as_array(value_at_pos_taag_targs)
_x_taat_wt_targs, _y_taat_wt_targs = as_array(value_at_pos_taat_wt_targs)
_x_taag_r1_targs, _y_taag_r1_targs = as_array(value_at_pos_taag_r1_targs)

matplotlib.rcParams.update({'font.size': 8})
ax = plt.subplot(2,2,1)
plt.xlim([-500,200]); plt.ylim([1,6])
plt.scatter(_x_taat_targs, np.log10(_y_taat_targs), alpha=0.05, edgecolor='none', color='black')
ax = plt.subplot(2,2,2)
plt.xlim([-500,200]); plt.ylim([1,6])
plt.scatter(_x_taag_targs, np.log10(_y_taag_targs), alpha=0.05, edgecolor='none', color='black')
ax = plt.subplot(2,2,3)
plt.xlim([-500,200]); plt.ylim([1,6])
plt.scatter(_x_taat_wt_targs, np.log10(_y_taat_wt_targs), alpha=0.05, edgecolor='none', color='black')
ax = plt.subplot(2,2,4)
plt.xlim([-500,200]); plt.ylim([1,6])
plt.scatter(_x_taag_r1_targs, np.log10(_y_taag_r1_targs), alpha=0.05, edgecolor='none', color='black')
plt.show()
