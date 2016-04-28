"""

Build a figure showing the effect of different controls
on motif enrichment.

Uses the columns added by:
filter_peaks_by_ratio_in_peak_region.py --load_bed_info

Expect three bar charts, sharing a y axis to the left of both,
denoting the contols. The first bar chart shows motif enrichment,
the second peak number. The third is the fraction of known
positives recovered.

"""
import sys
sys.path.insert(0, '/Network/Servers/file.biochem.wisc.edu/Volumes/BioSAN/Users/dfporter/.local/lib/python2.6/site-packages')
sys.path.insert(0, './src/')
#import filter_peaks_by_ratio_in_peak_region
import HTSeq
import pandas
import os
import argparse
import subprocess
import callpeaks
from callpeaks import get_sequences, score_binding_site
import matplotlib.pyplot as plt
from numpy.random import rand
from numpy import arange
import numpy as np
import re


def filter_list_by_ratio(peaks, col1, col2, min_1_to_2_ratio):
    filtered_peaks = []
    #print peaks
    for n, row in enumerate(peaks):
        if not n % 100: print "\t\col1 {ca}: {cav} col2 {cb}: {cbv}".format(
            ca=col1, cav=row[col1], cb=col2, cbv=row[col2]
        )
        denom = max(1., float(row[col2]))
        ratio = float(row[col1])/denom
        if ratio >= min_1_to_2_ratio:
            filtered_peaks.append(row)
    return filtered_peaks


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


def get_sequences(combined):
    return_df=False
    if type(combined) == pandas.DataFrame:
        combined = convert_df_to_list_of_dicts(combined)
        return_df=True
    #fasta_filename = '/home/dp/Desktop/celegans_genome/wormbase_ws235/c_elegans.WS235.genomic.fa'
    fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for peak_row in combined:
        start = peak_row['left']
        end = peak_row['right']
        chrm = peak_row['chrm']
        seq = sequences[chrm][start:end]
        #print "%s:%i-%i: seq %s" % (chrm, start, end, seq)
        if peak_row['strand'] == '-':
            seq = rc(seq)
        peak_row['seq'] = seq
    return pandas.DataFrame(combined)


def score_binding_site(peaks):
    for row in peaks:
        #print 'searching %s' % peaks.loc[index, 'seq']
        if re.search('tgt\w\w\wat', row['seq'].lower()) is not None:
            row['has_fbe'] = 1
            row['number_of_fbes'] = len(
                re.findall('tgt\w\w\wat', row['seq']))
        else:
            row['has_fbe'] = 0
            row['number_of_fbes'] = 0
        if re.search('ctgt\w\w\wat', row['seq'].lower()) is not None:
            row['minus_one_c'] = 1
        else:
            row['minus_one_c'] = 0
        if re.search('c\wtgt\w\w\wat', row['seq'].lower()) is not None:
            row['minus_two_c'] = 1
        else:
            row['minus_two_c'] = 0
        if row['minus_one_c'] == 1 or row['minus_two_c'] == 1:
            row['minus_one_or_two_c'] = 1
        else:
            row['minus_one_or_two_c'] = 0


def get_score(peaks_l):
    if type(peaks_l) == type(pandas.DataFrame()):
        peaks = peaks_l
        if ('seq' in peaks.columns) and ('has_fbe' in peaks.columns):
            pos = sum([1. for x in peaks.has_fbe if x>0])
            perc = 100. * float(pos)/float(len(peaks.index))
            return float(len(peaks.index)), perc
        peaks_l = convert_df_to_list_of_dicts(peaks_l)
    pos = neg = 0.
    if len(peaks_l) == 0: return (0, 0)
    if 'seq' not in peaks_l[0]:
        get_sequences(peaks_l)
    if 'has_fbe' not in peaks_l[0]:
        score_binding_site(peaks_l)
    for row in peaks_l:
        if row['has_fbe'] > 0: pos += 1.
        else: neg += 1.
    if int(pos + neg) > 0:
        perc = 100. * pos/(pos + neg)
    else: perc = 0
    num_peaks = pos + neg
    return (num_peaks, perc)


def convert_df_to_list_of_dicts(peaks):
    peaks_l = []
    for row in peaks.itertuples():
        peaks_l.append(
            dict((str(peaks.columns[num-1]), val) for num, val in enumerate(row))
            )
    return peaks_l


def score_binding_site_vs_cutoff(peaks, output_dirname='figs/',
                                 neg_control='n2'):
    li = ''
    results = []
    if type(peaks) == type(pandas.DataFrame()):
        peaks_l = convert_df_to_list_of_dicts(peaks)
    for base in range(0, 105, 5):
        base = float(base)
        ratio1 = base
        ratio2 = base# + 10
        if 'ratio' in peaks.columns:
            above = peaks[peaks['ratio']>=base]
            (num_peaks, perc_fbe) = get_score(above)
            perc_known_pos = score_positives(above)
        else:
            (num_peaks, perc_fbe, perc_known_pos) = filter_and_get_score(
                peaks_l, ratio1, ratio2, neg_control=neg_control)
        li += '%i peaks, %.1f percent with FBE.\n' % (int(num_peaks), perc_fbe)
        results.append((ratio1, ratio2, num_peaks, perc_fbe, perc_known_pos))
    print li
    y_tick_labels = tuple([str(x[0]) for x in results])
    plot_barchart(results, y_tick_labels, output_dirname=output_dirname, neg_control=neg_control)
    plot_results_line_graph(results, output_dirname=output_dirname, neg_control=neg_control)


def plot_barchart(results, y_tick_labels, output_dirname='figs/',
                  neg_control='unknown'):
    if not os.path.exists(output_dirname):
        os.system(output_dirname)
    plt.clf()
    pos = np.arange(len(results))
    print 'size of pos: %i' % len(pos)
    print 'pos: %s' % str(pos)
    print 'ratios: %s' % str([x[0] for x in results])
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
    #f, ax1 = plt.subplots(1, 1, sharey=True)
    # Plot 1.
    width = 0.8
    ax1.set_ylabel('Minimum FBF enrichment(reads/peak region ratio)',
                   fontsize=8)
    ax1.barh(pos, [x[2] for x in results], width, color='k')
    ax1.set_yticks(pos + 0.4)
    ax1.set_yticklabels(y_tick_labels)
    max_peak_num = sorted([x[2] for x in results])[-1]
    print sorted([x[2] for x in results])
    print max_peak_num
    print 'ticks put at %s' % str(range(0, int(max_peak_num) + 1000, 1000))
    if max_peak_num > 6000.:
        ax1.set_xticks(
            range(0., max_peak_num + 2000., 2000.))
    else:
        ax1.set_xticks(
            range(0., max_peak_num + 1000., 1000.))
    ax1.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='y', labelsize=8)
    ax1.axvline(0, color='k', lw=3)
    ax1.set_xlabel('Number of peaks', fontsize=8)
    ax1.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off')
    ax1.xaxis.grid(True)
    #plt.plot()
    ax1.set_ylim([0, len(results)])
    plt.savefig('{subdir}/effect_of_controls_barplot_{nc}.pdf'.format(
        subdir=output_dirname, nc=neg_control), format='pdf')
    #sys.exit()
    # Plot 2.
    ax2.set_yticks(pos, tuple([str(x[0]) for x in results]))
    ax2.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off')
    ax2.barh(pos, [x[3] for x in results], width, color='k')
    ax2.axvline(0, color='k', lw=3)
    ax2.set_xlabel('Percent of peaks with FBE', fontsize=8)
#    plt.setp(ax2.get_xticklabels(), fontsize=7)
    ax2.tick_params(axis='x', labelsize=8)
    ax2.xaxis.grid(True)
    # Plot 3.
    ax3.set_yticks(pos, tuple([str(x[0]) for x in results]))
    ax3.barh(pos, [x[4] for x in results], width, color='k')
    ax3.set_xticks(range(0, 120, 20))
    ax3.axvline(0, color='k', lw=3)
    ax3.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off')
    ax3.set_xlabel('Percent of known positives identified', fontsize=8)
    ax3.tick_params(axis='x', labelsize=8)
    ax3.xaxis.grid(True)
#    plt.setp(ax3.get_xticklabels(), fontsize=7)
    plt.tight_layout(pad=1)
    plt.savefig('{subdir}/effect_of_controls_barplot_{nc}.pdf'.format(
        subdir=output_dirname, nc=neg_control), format='pdf')
    plt.clf()


def plot_results_line_graph(results, output_dirname='figs/',
                            neg_control='unknown'):
    if not os.path.exists(output_dirname):
        os.system(output_dirname)
    plt.clf()
    plt.xlabel('Minimum FBF enrichment (reads/peak region ratio).')
    plt.ylabel('Number of peaks.')
    plt.plot([x[0] for x in results],
             [x[2] for x in results],
             label='Number of peaks')
    plt.savefig('{subdir}/min_ratio_vs_num_peaks.pdf'.format(
        subdir=output_dirname), format='pdf')
    plt.clf()
    plt.xlabel('Minimum FBF enrichment (reads/peak region ratio).')
    plt.ylabel('Peaks with FBE (\%).')
    plt.plot([x[0] for x in results],
             [x[3] for x in results],
             label='Peaks with FBE')
    plt.savefig('{subdir}/min_ratio_vs_perc_fbe.pdf'.format(
        subdir=output_dirname), format='pdf')
    plt.clf()
    plt.xlabel('Number of peaks.')
    plt.ylabel('Peaks with FBE (\%)')
    plt.plot([x[2] for x in results],
             [x[3] for x in results],
             label='Number of peaks vs presence binding site .')
    plt.savefig('{subdir}/perc_fbe_vs_num_peaks_{nc}.pdf'.format(
        subdir=output_dirname, nc=neg_control), format='pdf')
    plt.clf()


def filter_and_get_score(peaks_l, ratio1, ratio2, neg_control='n2',
                         normalize_fbf1=1, normalize_fbf2=1, only_use=['fbf1']):
    if neg_control=='n2':
            if 'fbf1' in only_use:
                peaks_l = filter_list_by_ratio(
                    peaks_l, 'fbf1_reads', 'n2_from_fbf1_reads', ratio1)
            if 'fbf2' in only_use:
                peaks_l = filter_list_by_ratio(
                    peaks_l, 'fbf2_reads', 'n2_from_fbf2_reads', ratio2)
    if neg_control=='rna_seq':
            if 'fbf1' in only_use:
                peaks_l = filter_list_by_ratio(
                    peaks_l, 'fbf1_reads', 'rna_seq_reads', float(ratio1)/float(normalize_fbf1))
            if 'fbf2' in only_use:
                peaks_l = filter_list_by_ratio(
                    peaks_l, 'fbf2_reads', 'rna_seq_reads', float(ratio2)/float(normalize_fbf2))
    (perc, num_peaks) = get_score(peaks_l)
    perc_knowns = score_positives(peaks_l)
    return (perc, num_peaks, perc_knowns)

def score_positives(peaks_l):
    known_pos = set(['gld-1', 'htp-1', 'htp-2', 'mpk-1', 'him-3', 'fbf-1', 'lip-1', 'syp-2', 'fbf-2', 'fog-1', 'fem-3', 'syp-3', 'gld-3', 'fog-3', 'egl-4'])
    if type(peaks_l) == type([]):
        obs_genes = set([x['gene_name'] for x in peaks_l])
    if type(peaks_l) == type(pandas.DataFrame(peaks_l)):
        obs_genes = set(peaks_l['gene_name'].tolist())
    obs_pos = known_pos & obs_genes
    missing_pos = known_pos - obs_genes
    print "\tObserved %s.\n\tMissing %s..." % (str(obs_pos), str(missing_pos))
    obs_pos_n = len(list(obs_pos))
    missing_pos_n = len(list(missing_pos))
    perc_knowns = 100. * float(obs_pos_n)/float(len(list(known_pos)))
    return perc_knowns


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Make a figure of the effect of different controls
on FBE enrichment and peak numbers.
Uses the columns added by:
filter_peaks_by_ratio_in_peak_region.py --load_bed_info.
""")
    parser.add_argument('-i', '--input', dest='input',
                        default='with_peak_nums/combined_fbf_4_rna_seq.txt',
                        help='Peaks filename.')
    args = parser.parse_args()
    peaks_filename = args.input
    peaks = pandas.read_csv(peaks_filename, sep='\t', index_col=False)
    print peaks
#    out_filename = 'with_peak_nums/'
#    out_filename += os.path.basename(peaks_filename).partition('.txt')[0]
    # Bar 1: vs. neg IP.
#    peaks = convert_df_to_list_of_dicts(peaks)
#    get_sequences(peaks)
#    score_binding_site(peaks)
    score_binding_site_vs_cutoff(
        peaks, output_dirname='figs/methods_comparison')

