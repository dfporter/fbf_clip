import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
#import pysam
import sys
import argparse
import matplotlib
import os

def make_fig(peaks, output_filename='boxplot_of_height_v_position.pdf', split_y=None,
             max_y=1e3):
    if split_y is not None:
        # Good value is 1500.
        with_split(peaks, output_filename=output_filename, split_y=split_y,
                   max_y=max_y)
        return
    vals = {}
    for index, row in peaks.iterrows():
        vals.setdefault(row['location'], [])
        vals[row['location']].append(row['height'])
    labels = ['ncRNA',  """5'UTR""",'CDS', """3'UTR""", ]
    labels = [label for label in labels if label in vals]
    list_of_vals = [vals[label] for label in labels]
    fig, ax1 = plt.subplots(1, 1, sharex=True)
    ax1.set_ylim(0, max_y)
    bp1 = ax1.boxplot(list_of_vals,labels=labels)
    plt.setp(bp1['boxes'], color='black', linewidth=2)
    plt.setp(bp1['whiskers'], color='black', linewidth=2)
    plt.setp(bp1['fliers'], color='black', linewidth=1, marker='o', alpha=0.2)
    plt.setp(bp1['medians'], color='black', linewidth=2)
#    ax2.set_ylim([0, split_y])
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.xaxis.tick_bottom()

    pos = np.arange(len(labels)) + 1
    for tick, label in zip(range(len(labels)), labels):
        print tick
        print label
        k = tick % 2
        ax1.text(
            pos[tick], max_y*0.90, 'n=' + str(len(vals[label])),
            horizontalalignment='center')
    fig.text(0.04, 0.5, 'Peak height', va='center', rotation='vertical')
    plt.savefig(output_filename, format='pdf')
    plt.clf()

def with_split(peaks, output_filename='boxplot_of_height_v_position.pdf',
               split_y=1500, max_y=1000):
    vals = {}
    high_vals = {}
    low_vals = {}
    for index, row in peaks.iterrows():
        high_vals.setdefault(row['location'], [])
        low_vals.setdefault(row['location'], [])
        vals.setdefault(row['location'], [])
        vals[row['location']].append(row['height'])
        if row['height'] > 100:
            high_vals[row['location']].append(row['height'])
        else:
            low_vals[row['location']].append(row['height'])
    labels = ['ncRNA', """3'UTR""", """5'UTR""", 'CDS']
    labels = [label for label in labels if label in vals]
    list_of_vals = [vals[label] for label in labels]
    list_of_high_vals = [high_vals[label] for label in labels]
    list_of_low_vals = [low_vals[label] for label in labels]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    bp1 = ax1.boxplot(list_of_vals,labels=labels)
    bp2 = ax2.boxplot(list_of_vals,labels=labels)
    plt.setp(bp1['boxes'], color='black', linewidth=2)
    plt.setp(bp2['boxes'], color='black', linewidth=2)
    plt.setp(bp1['whiskers'], color='black', linewidth=2)
    plt.setp(bp2['whiskers'], color='black', linewidth=2)
    plt.setp(bp1['fliers'], color='black', linewidth=1, marker='o', alpha=0.2)
    plt.setp(bp2['fliers'], color='black', linewidth=1, marker='o', alpha=0.2)
    plt.setp(bp1['medians'], color='black', linewidth=2)
    plt.setp(bp2['medians'], color='black', linewidth=2)
    ax1.set_ylim([split_y, max_y])
    ax2.set_ylim([0, split_y])
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop='off')
    ax2.xaxis.tick_bottom()

    d = .005 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d,+d),(-d,+d), **kwargs)      # top-left diagonal
    ax1.plot((1-d,1+d),(-d,+d), **kwargs)    # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d),(1-d,1+d), **kwargs)   # bottom-left diagonal
    ax2.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-right diagonal

    pos = np.arange(len(labels)) + 1
    top = 1000
    for tick, label in zip(range(len(labels)), labels):
        print tick
        print label
        k = tick % 2
        ax1.text(
            pos[tick], top*0.90, 'n=' + \
            str(len(low_vals[label]) + len(high_vals[label])),
            horizontalalignment='center')
    fig.text(0.04, 0.5, 'Peak height', va='center', rotation='vertical')
    plt.savefig(output_filename, format='pdf')
    plt.clf()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input',
        default='/groups/Kimble/Aman Prasad/clip/methods/filtered_gauss_nil_01_five_reps/combined_fbf.txt',
        help='''Input file of peaks, 5/6 reps.''')
    parser.add_argument(
        '-s', '--split_y',
        default=None,
        help='''Split y axis.''')
    args = parser.parse_args()
    max_y = 10000.
    if args.split_y is not None:
        args.split_y = float(args.split_y)
        if args.split_y >= (max_y * 0.8):
            args.split_y = max_y/2.
    #else:
    #    args.split_y = max_y/2.
    #filename = 'with_info_pypeaks_fdr1_negip_local/combined_fbf1.txt'
    peaks = pandas.read_csv(args.input, sep='\t')
    if ('height' not in peaks.columns) or (
        len(peaks['height'].value_counts()) == 1):
        print "%s seems to be lacking a height column with more \
than one value." % args.input
        if raw_input(
            "Add height as fbf1_reads + fbf2_reads? \
This will be a sum of the individual reads per million values. (Y/N):").upper() == 'Y':
            tups = zip(peaks['fbf1_reads'].tolist(), peaks['fbf2_reads'].tolist())
            peaks['height'] = [sum(list(x)) for x in tups]
            peaks.to_csv(args.input, sep='\t', index=False)
    out_dir = os.path.dirname(args.input)
    if not os.path.exists('figs/' + out_dir):
        os.system('mkdir figs/' + out_dir)
    output_filename = 'figs/{d}/{f}'.format(
        d=out_dir,
        f='boxplot_of_height_v_position.pdf')
    print output_filename
    make_fig(peaks, output_filename=output_filename, split_y=args.split_y,
             max_y=max_y)
