"""

1. Score position of all FBE binding sites in each sequence.
2. Score the position of the deletion (assume center nucleotide).
3. Add a column of deletion - FBE start.
4. Create a histogram by position: freq = f(dist)
5. plot the histogram as a line.

"""


import os
import sys
import re
import pandas
import glob
from get_sequences_for_table import get_sequences_for_table
import HTSeq
import matplotlib
import matplotlib.pyplot as plt
import operator


def cims_pos_vs_fbe(filename1, filename2, sequences,
                    two_freq_filename='../clip/figs/cim_vs_fbe_separate.pdf',
                    combined_filename='../clip/figs/cim_vs_fbe_fbf1_and_2.pdf'
                    ):
    hist1 = score_dist_to_fbe(filename1, sequences)#, pat_seq='ACA\w\w\wTA')
    hist2 = score_dist_to_fbe(filename2, sequences)#, pat_seq='ACA\w\w\wTA')
    combined = {}
    for dist in hist1:
        hist2.setdefault(dist, 0)
        combined[dist] = hist1[dist] + hist2[dist]
    for dist in hist2:
        if dist in hist1: continue
        combined[dist] = hist2[dist]
    print "Combined: %s" % str(combined)
    combined = sorted(combined.items(), key=operator.itemgetter(0))
    hist1 = sorted(hist1.items(), key=operator.itemgetter(0))
    hist2 = sorted(hist2.items(), key=operator.itemgetter(0))
    total_y = sum([t[1] for t in combined])
    total_1 = sum([t[1] for t in hist1])
    total_2 = sum([t[1] for t in hist2])
    print "Total instances: %i" % total_y
    print "Total instances, FBF-1: %i" % total_1
    print "Total instances, FBF-2: %i" % total_2
#    plot_freq(hist1, '../clip/figs/cim_vs_fbe_fbf1.pdf')
#    plot_freq(hist2, '../clip/figs/cim_vs_fbe_fbf2.pdf')
    plot_two_freq(hist1, hist2, two_freq_filename)
    plot_freq(combined, combined_filename)


def score_dist_to_fbe(filename, sequences, pat_seq='TGT\w\w\wAT'):
    verbose = False
    peaks = pandas.read_csv(filename, sep='\t')
#    if 'seq' not in peaks.columns:
    get_sequences_for_table(
        peaks, sequences, expand=20)
    pat = re.compile(pat_seq)
    motif_pos = {}
    if len(peaks['seq']) < 2:
        return
    hist = {}
    nt_count = {}
    for index, row in peaks.iterrows():
        center_nt = len(peaks.loc[index, 'seq'])/2
        motif_pos[index] = []
        for m in pat.finditer(peaks.loc[index, 'seq']):
            dist = center_nt - m.start()
            nt_count.setdefault(dist,
                                {'A':0, 'T': 0, 'G': 0, 'C':0})
            xl_nt = peaks.loc[index, 'seq'][center_nt:center_nt+1]
            nt_count[dist][
                peaks.loc[index, 'seq'][center_nt:center_nt+1]
                ] += 1
            if abs(dist) < 2:
                if verbose: print "xl dist %i: %s" % (
                    dist,
                    peaks.loc[index, 'seq'][center_nt:center_nt+2])
                if verbose: print "\tFBE seq: %s" % (
                    peaks.loc[index, 'seq'][m.start():m.start()+8])
                if verbose: print "\tSeq: %s" % peaks.loc[index, 'seq']
                if verbose: print "\tcenter_nt: %i m.start(): %i dist: %i" % (
                    center_nt, m.start(), dist)
            motif_pos[index].append({
                'pos': m.start(),
                'dist': dist})
            hist.setdefault(dist, 0)
            hist[dist] += 1
    for pos in nt_count:
        if verbose: print "%i: %s" % (pos, str(nt_count[pos]))
    meme_filename = '../clip/figs/'
    meme_filename += os.path.basename(filename).partition('.')[0] + '.meme'
    #print_seqlogo_format(nt_count, out_filename=meme_filename)
    print_meme_format(nt_count, out_filename=meme_filename)
    if verbose: print hist
    return hist

def print_meme_format(nt_count, out_filename='../clip/figs/logo.meme'):
    start = -10
    end = 11
    outli = """MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF cims
letter-probability matrix: alength= 4 w= %i nsites= 100 E=1e-10
""" % int(end - start -1)
    for pos in range(start, end):
        if pos in nt_count:
            total = float(sum([nt_count[pos][x] for x in nt_count[pos]]))
            outli += """%f\t%f\t%f\t%f\n""" % (
                float(nt_count[pos]['A'])/total,
                float(nt_count[pos]['C'])/total,
                float(nt_count[pos]['G'])/total,
                float(nt_count[pos]['T'])/total)
        else:
            outli += """%f\t%f\t%f\t%f\n""" % (.25,.25,.25,.25)
    with open(out_filename, 'w') as f:
        f.write(outli)
    logo_filename = '../clip/figs/%s' % os.path.basename(
        out_filename).partition('.')[0] + '_pwm_logo.eps'
    os.system('ceqlogo -i %s -o %s -f %s' % (out_filename, logo_filename,
                                       'EPS'))
    tmp_file = open('tmp.eps', 'w')
    with open(logo_filename, 'r') as f:
        for li in f:
            li = re.sub('\(T\) numchar', '(U) numchar', li)
            tmp_file.write(li)
    os.system('mv tmp.eps' + ' ' + logo_filename)


def print_seqlogo_format(nt_count, out_filename='../clip/figs/logo.seqlogo'):
    start = -10
    end = 11
    outli = ""
    # 'A' on line 1.
    for base in ['A', 'C', 'G', 'T']:
        if base != 'T':
            outli += "%s |" % base
        else:  # T>U
            outli += "U |"
        for pos in range(start, end):
            if pos in nt_count:
                total = float(sum([nt_count[pos][x] for x in nt_count[pos]]))
                outli += """ %f""" % (
                    float(nt_count[pos][base])/total)
            else: outli += ' 0.25 0.25 0.25 0.25'
        outli += '\n'
    with open(out_filename, 'w') as f:
        f.write(outli)
    logo_filename = '../clip/figs/%s' % os.path.basename(
        out_filename).partition('.')[0] + '_pwm_seqlogo.pdf'
    rcmd = '''
library(seqLogo)
library(Biostrings)
pcm = read.table("%s")
pcm = pcm[,3:ncol(pcm)]
rownames(pcm) = c("A", "C", "G", "U")
pwm = makePWM(pcm) #alphabet='AA')
pdf("%s", width=10, height=5)
seqLogo(pcm, ic.scale=FALSE)# xfontsize=10, yfontsize=15)
dev.off()
''' % (out_filename, logo_filename)
    with open('tmp.r', 'w') as f:
        f.write(rcmd)
    os.system('R CMD BATCH tmp.r')
    #os.system('ceqlogo -i %s -o %s -f %s' % (out_filename, logo_filename,
    #                                   'EPS'))


def plot_freq(freq, output_filename, motif='UGUNNNAU'):
    plt.clf()
    total_y = sum([t[1] for t in freq])
    if total_y != 1:
        y = [100 * float(t[1])/float(total_y) for t in freq]
    else:
        y = [t[1] for t in freq]
    x = [t[0] for t in freq]
    plt.xlabel('Position relative to FBE (nt)')#, labelpad=20)
    plt.ylabel('Frequency (%)')
    plt.plot(x,y)
    for index, base in enumerate(list(motif)):
        plt.text(index,0.07, base, fontsize=15,
                 horizontalalignment='center',)
    plt.savefig(output_filename, format='pdf')
    plt.clf()

def plot_two_freq(freq1, freq2, output_filename, motif='UGUNNNAU'):
    plt.clf()
    total_y1 = sum([t[1] for t in freq1])
    if total_y1 != 1:
        y1 = [100 * float(t[1])/float(total_y1) for t in freq1]
    else:
        y1 = [t[1] for t in freq1]
    x1 = [t[0] for t in freq1]
    total_y2 = sum([t[1] for t in freq2])
    if total_y1 != 1:
        y2 = [100 * float(t[1])/float(total_y2) for t in freq2]
    else:
        y2 = [t[1] for t in freq2]
    x2 = [t[0] for t in freq2]
    plt.xlabel('Position relative to FBE (nt)')
    plt.ylabel('Frequency (%)')

    plt.plot(x1,y1, 'b-', label='FBF-1')
    plt.plot(x2,y2, 'r-', label='FBF-2')
    for index, base in enumerate(list(motif)):
        plt.text(index,0.07, base, fontsize=15,
                 horizontalalignment='center')
    plt.legend()
    plt.savefig(output_filename, format='pdf')
    plt.clf()

def cim_vs_site_in_dir(in_dir):
    fasta_filename = '/scratch/indexes/WS235.fa'
    sequences = dict(
        (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
    #gtf_df = pandas.read_csv('../clip/lib/gtf_with_names_column.txt', sep='\t')
    #sk = '''
    cims_pos_vs_fbe(
        '%s/combined_fbf1_cits.txt' % in_dir,
        '%s/combined_fbf2_cits.txt' % in_dir,
        sequences,
        two_freq_filename='../clip/figs/cits_pos_vs_fbe.pdf',
        combined_filename='../clip/figs/cits_combined_pos_vs_fbe.pdf')
    cims_pos_vs_fbe(
        '%s/combined_fbf1_cims.txt' % in_dir,
        '%s/combined_fbf2_cims.txt' % in_dir,
        sequences,
        two_freq_filename='../clip/figs/cims_pos_vs_fbe.pdf',
        combined_filename='../clip/figs/cims_combined_pos_vs_fbe.pdf')
    #'''
    for filename in [
        '%s/combined_fbf1_cits.txt' % in_dir,
        '%s/combined_fbf2_cits.txt' % in_dir,
        '%s/combined_fbf1_cims.txt' % in_dir,
        '%s/combined_fbf2_cims.txt' % in_dir,]:
        score_dist_to_fbe(filename, sequences)#, pat_seq='ACA\w\w\wTA')


if __name__ == '__main__':
    cim_vs_site_in_dir('reproducible/')
