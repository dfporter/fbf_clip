import glob
import os
import re
import sys
import matplotlib.pyplot as plt
import pandas
from get_sequences_for_table import get_sequences_for_table
import HTSeq
import argparse

def add_minus_three_c_column(peaks):
    if 'seq' not in peaks.columns:
        fasta_filename = '/scratch/indexes/WS235.fa'
        sequences = dict(
            (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
        get_sequences_for_table(peaks, sequences, expand=10)
    peaks['minus_three_c'] = 0
    peaks['minus_four_c'] = 0
    peaks['tgt'] = 0
    peaks['has_fbe'] = 0
    peaks['seq'] = [x.lower() for x in peaks['seq'].tolist()]
    for index, row in peaks.iterrows():
        if re.search('tgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'has_fbe'] = 1
        else:
            peaks.loc[index, 'has_fbe'] = 0
        if re.search('c\w\wtgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'minus_three_c'] = 1
        else: peaks.loc[index, 'minus_three_c'] = 0
        if re.search('c\w\w\wtgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'minus_four_c'] = 1
        else: peaks.loc[index, 'minus_four_c'] = 0
        if re.search('tgt', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'tgt'] = 1
        else: peaks.loc[index, 'tgt'] = 0
        if re.search('ctgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'minus_one_c'] = 1
        else: peaks.loc[index, 'minus_one_c'] = 0
        if re.search('c\wtgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'minus_two_c'] = 1
        else: peaks.loc[index, 'minus_two_c'] = 0
    return peaks

def add_cumulative_column(peaks, use_col='height'):
    print "Currently broken - need height"
    peaks.sort(use_col, inplace=True, ascending=False)
    peaks['frac_fbe'] = 0
    peaks['frac_minus_one_c'] = 0
    peaks['frac_minus_two_c'] = 0
    peaks['frac_minus_three_c'] = 0
    peaks['frac_minus_four_c'] = 0
    peaks['frac_tgt'] = 0
    peaks['rank'] = 1
    prev_rank = 0
    sum_fbe = 0.0
    sum_with_minus_one = 0.0
    sum_with_minus_two = 0.0
    sum_with_minus_three = 0.0
    sum_with_minus_four = 0.0
    sum_with_tgt = 0
    sum_rows = 0.0
    for index, row in peaks.iterrows():
        sum_rows += 1.0
        peaks.loc[index, 'rank'] = prev_rank + 1
        prev_rank += 1
        #print "Index %i, height %f" % (
        #    index, float(peaks.loc[index, 'height']))
        if peaks.loc[index, 'has_fbe'] > 0:
            sum_fbe += 1.0
        if peaks.loc[index, 'minus_one_c'] > 0:
            sum_with_minus_one += 1.0
        if peaks.loc[index, 'minus_two_c'] > 0:
            sum_with_minus_two += 1.0
        if peaks.loc[index, 'minus_three_c'] > 0:
            sum_with_minus_three += 1.0
        if peaks.loc[index, 'minus_four_c'] > 0:
            sum_with_minus_four += 1.0
        if peaks.loc[index, 'frac_tgt'] > 0:
            sum_with_tgt += 1.0
        peaks.loc[index, 'frac_fbe'] = sum_fbe/sum_rows
        peaks.loc[index, 'frac_minus_one_c'] = sum_with_minus_one/sum_rows
        peaks.loc[index, 'frac_minus_two_c'] = sum_with_minus_two/sum_rows
        peaks.loc[index, 'frac_minus_three_c'] = sum_with_minus_three/sum_rows
        peaks.loc[index, 'frac_minus_four_c'] = sum_with_minus_four/sum_rows
        peaks.loc[index, 'frac_tgt'] = sum_with_tgt/sum_rows
    return peaks

def save_new_file(peaks, folder_name):
    if not os.path.exists(folder_name):
        os.system('mkdir ' + folder_name)


def load_and_annotate_peaks(fname, use_col='', norm_to_rna_seq=False,
                            norm_col_name=''):
    peaks = pandas.read_csv(fname, sep='\t')
    peaks = peaks[peaks['biotype']=='protein_coding']
    peaks = add_minus_three_c_column(peaks)
    print fname, use_col, norm_to_rna_seq
    print peaks.head(1)
    print peaks.tail(1)
    if norm_to_rna_seq:
        min_rna_seq = min([x for x in peaks['rna_seq_reads'].tolist() if x>0])
        tups = zip(peaks[use_col].tolist(), peaks['rna_seq_reads'].tolist())
        tups = [(x[0], max([float(x[1]), min_rna_seq])) for x in tups]
        print tups
        peaks[norm_col_name] = [
            float(x[0])/float(x[1]) for x in tups]
        peaks = add_cumulative_column(peaks, use_col=norm_col_name)
        return peaks
    return add_cumulative_column(peaks, use_col=use_col)

    
def plot_fbes_vs_rank(peaks1, peaks2, out_filename, combined=True):
    if not combined:
        plt.clf()
        plt.plot(peaks1['rank'], peaks1['frac_fbe'], 'b-',
                 linewidth=3, label='FBF-1 with FBE')
        plt.plot(peaks1['rank'], peaks1['frac_minus_one_c'], 'g-',
                 linewidth=3,  label='FBF-1 -1C')
        plt.plot(peaks1['rank'], peaks1['frac_minus_two_c'], 'm-',
                 linewidth=3,  label='FBF-1 -2C')
        plt.plot(peaks1['rank'], peaks1['frac_minus_three_c'], 'c-', label='-3C')
        plt.plot(peaks1['rank'], peaks1['frac_minus_four_c'], 'k-', label='-4C')
        plt.legend(prop={'size':12, 'weight': 'bold'})
        plt.xlabel('Peak rank by height')
        plt.ylabel('Fraction of peaks with given motif')
        plt.savefig(out_filename)
        plt.clf()
        return
    plt.clf()
    plt.plot(peaks1['rank'], peaks1['frac_fbe'], 'b-',
             linewidth=3, label='FBF-1 with FBE')
    plt.plot(peaks1['rank'], peaks1['frac_minus_one_c'], 'g--',
             linewidth=3,  label='FBF-1 -1C')
    plt.plot(peaks2['rank'], peaks2['frac_fbe'], 'r-',
             linewidth=3, label='FBF-2 with FBE')
    plt.plot(peaks2['rank'], peaks2['frac_minus_one_c'], 'c--',
             linewidth=3,  label='FBF-2 -1C')
    plt.plot(peaks1['rank'], peaks1['frac_minus_two_c'], 'm--',
             linewidth=3,  label='FBF-1 -2C')
    plt.plot(peaks2['rank'], peaks2['frac_minus_two_c'], 'k--',
             linewidth=3, label='FBF-2 -2C')
#    plt.plot(peaks['rank'], peaks['frac_minus_three_c'], 'r-', label='-3C')
#    plt.plot(peaks['rank'], peaks['frac_minus_four_c'], 'm-', label='-4C')
    plt.legend(prop={'size':12, 'weight': 'bold'})
    plt.xlabel('Peak rank by height')
    plt.ylabel('Fraction of peaks with given motif')
    if not os.path.exists(os.path.dirname(out_filename)):
        os.system('mkdir ' + os.path.dirname(out_filename))
    plt.savefig(out_filename)
    plt.clf()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default='peaks/',
        help='Directory folder holding combined_fbf*.txt peaks files.')
    parser.add_argument('-n', '--norm_to_rna_seq',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    input_dirname = args.input
    label = input_dirname.rstrip('/') #os.path.dirname(input_filename)
    if label == '.' or label == '':
        label = 'no_label'
    if not os.path.exists('figs/' + label + '/'):
        os.system('mkdir figs/' + label)    
    if args.norm_to_rna_seq:
        peaks1 = load_and_annotate_peaks(input_dirname + '/combined_fbf1.txt',
                               use_col='fbf1_reads',
                               norm_col_name='fbf1_reads_over_rna_seq',
                               norm_to_rna_seq=args.norm_to_rna_seq)
        peaks2 = load_and_annotate_peaks(input_dirname + '/combined_fbf2.txt',
                               use_col='fbf2_reads',
                               norm_col_name='fbf2_reads_over_rna_seq',
                               norm_to_rna_seq=args.norm_to_rna_seq)
        plot_fbes_vs_rank(
            peaks1, peaks1, 'figs/{lab}/{base}_one_to_four.eps'.format(
            lab=label, base='fbf1' + '_norm_to_rna_seq'), combined=False)
        plot_fbes_vs_rank(
            peaks2, peaks2, 'figs/{lab}/{base}_one_to_four.eps'.format(
            lab=label, base='fbf2' + '_norm_to_rna_seq'), combined=False)
        plot_fbes_vs_rank(
            peaks1, peaks2, 'figs/{lab}/{base}_one_to_two.eps'.format(
            lab=label, base='separate' + '_norm_to_rna_seq'))
    peaks1 = load_and_annotate_peaks(input_dirname + '/combined_fbf1.txt',
                           use_col='fbf1_reads',
                           norm_to_rna_seq=False)
    peaks2 = load_and_annotate_peaks(input_dirname + '/combined_fbf2.txt',
                           use_col='fbf2_reads',
                           norm_to_rna_seq=False)
    plot_fbes_vs_rank(peaks1, peaks1, 'figs/{lab}/{base}_one_to_four.eps'.format(
        lab=label, base='fbf1'), combined=False)
    plot_fbes_vs_rank(peaks2, peaks2, 'figs/{lab}/{base}_one_to_four.eps'.format(
        lab=label, base='fbf2'), combined=False)
    plot_fbes_vs_rank(peaks1, peaks2, 'figs/{lab}/{base}_one_to_two.eps'.format(
        lab=label, base='separate'))
