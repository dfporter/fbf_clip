import pandas
import csv
import re
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
#import pysam
import matplotlib
import glob
import sys
import os
from scipy import stats
import argparse
import re


def rename_columns(peaks, which='both'):
    if which == 'fbf1':
        peaks['clip_reads'] = peaks['fbf1_reads']
        peaks['neg_ip_reads'] = peaks['n2_from_fbf1_reads']
    elif which == 'fbf2':
        peaks['clip_reads'] = peaks['fbf2_reads']
        peaks['neg_ip_reads'] = peaks['n2_from_fbf2_reads']
    elif which == 'both':
        peaks['clip_reads'] = float(peaks['fbf1_reads']) + float(peaks['fbf2_reads'])
        peaks['neg_ip_reads'] = float(peaks['n2_from_fbf1_reads']) + float(peaks['n2_from_fbf2_reads'])
    else:
        print "rename_columns 'which' option not set correctly (value=%s)" % str(which)
    tup = zip(peaks['clip_reads'].tolist(), peaks['neg_ip_reads'].tolist())
    peaks['clip_over_neg_ip'] = [float(x[1])/float(x[0]) for x in tup]
    return peaks


def plot_tup_heatmap_on_axis(tups, ax, bins=50):
    tups = [x for x in tups if (np.isfinite(np.log10(x[0]))) and (
        np.isfinite(np.log10(x[1])))]
    heatmap, xedges, yedges = np.histogram2d(
        [np.log10(x[0]) for x in tups], [np.log10(x[1]) for x in tups],
        bins=bins)
    min_low =  min([xedges[0], yedges[0]])  + 1
    max_high = max([xedges[-1], yedges[-1]]) - 1
    extent =  [min_low, max_high, min_low, max_high]
    extent = [min([x[0] for x in tups]), max([x[0] for x in tups]),
            min([x[1] for x in tups]), max([x[1] for x in tups])]
    extent = [np.log10(x) for x in extent]
    extent[0] = min([extent[0], extent[2]])
    extent[1] = max([extent[1], extent[3]])
    extent[2] = min([extent[0], extent[2]])
    extent[3] = max([extent[1], extent[3]])

    ax.hexbin([np.log10(x[0]) for x in tups],
              [np.log10(x[1]) for x in tups], extent=extent, bins='log')
    ax.axis(extent)
    return
    ax.autoscale(False)
    ax.axis(extent)
    ax.imshow(heatmap.T, extent=extent, origin='lower')
    #ax.set_xlim([-10, 10])
    #ax.set_ylim([-10, 10])
    
def scatterplot_vs_controls(peaks, label='', top_level_dir='', use_set_limits=True):
    top_level_dir = os.path.dirname(top_level_dir + '/')
    bins=50
    plt.clf()
    plt.figure(figsize=(8,8))
    matplotlib.rcParams.update({'font.size': 12})
    gs = matplotlib.gridspec.GridSpec(2,2)
    #plt.axes().set_aspect('equal')
    ax1 = plt.subplot(gs[0, 0], aspect='equal')
    ax1.set_title('Reads in peak regions')
    ax1.set_xlabel('$\mathregular{log_{10}}$ Negative IP reads')
    ax1.set_ylabel('$\mathregular{log_{10}}$ CLIP-seq reads')
    neg_max = np.round(np.max(np.log10(peaks.neg_ip_reads)), decimals=0)
    clip_max = np.round(np.max(np.log10(peaks.clip_reads)), decimals=0)
    rna_seq_max = np.max(np.log10(peaks.rna_seq_reads))
#    finite_ratios = peaks['clip_over_neg_ip']
#    finite_ratios = finite_ratios[finite_ratios<1e3]
#    ratio_max = np.max(np.log10(finite_ratios))
    tups = zip(peaks.neg_ip_reads, peaks.fbf1_reads)
    plot_tup_heatmap_on_axis(tups, ax1, bins=bins)
    (pearson_r, p_value) = stats.pearsonr(peaks.neg_ip_reads, peaks.clip_reads)
    ax1.text(0.6, 0.85, "$\mathregular{R^2}$ %.2f" % pearson_r**2, transform=ax1.transAxes)
    ax2 = plt.subplot(gs[0,1], aspect='equal')
    ax2.set_title('Reads in peak regions')
    ax2.set_ylabel('$\mathregular{log_{10}}$ CLIP-seq reads')
    ax2.set_xlabel('$\mathregular{log_{10}}$ RNA-seq reads')
    tups = zip(peaks.rna_seq_reads, peaks.clip_reads)
    plot_tup_heatmap_on_axis(tups, ax2, bins=bins)
    (pearson_r, p_value) = stats.pearsonr(peaks.rna_seq_reads, peaks.clip_reads)
    ax2.text(0.6, 0.85, "$\mathregular{R^2}$ %.2f" % pearson_r**2, transform=ax2.transAxes)
    ax3 = plt.subplot(gs[1,0], aspect='equal')
    ax3.set_title('Reads in peak regions')
    ax3.set_xlabel('$\mathregular{log_{10}}$ RNA-seq reads')
    ax3.set_ylabel('$\mathregular{log_{10}}$ Negative IP reads')
    tups = zip(peaks.rna_seq_reads, peaks.neg_ip_reads)
    plot_tup_heatmap_on_axis(tups, ax3, bins=bins)
    #ax3.set_xticks(np.arange(0, max([rna_seq_max, neg_max])+1.1, 1.0))
    #ax3.set_yticks(np.arange(0, max([rna_seq_max, neg_max])+1.1, 1.0))
    (pearson_r, p_value) = stats.pearsonr(peaks.neg_ip_reads, peaks.rna_seq_reads)
    ax3.text(0.6, 0.85, "$\mathregular{R^2}$ %.2f" % pearson_r**2, transform=ax3.transAxes)
    ax4 = plt.subplot(gs[1,1], aspect='equal')
    ax4.set_title('Reads in peak regions')
    ax4.set_xlabel('$\mathregular{log_{10}}$ RNA-seq reads')
    ax4.set_ylabel('$\mathregular{log_{10}}$ CLIP-seq/Negative IP')
    tups = zip(peaks.rna_seq_reads, peaks.clip_over_neg_ip)
    plot_tup_heatmap_on_axis(tups, ax4, bins=bins)
    #(pearson_r, p_value) = stats.pearsonr(peaks.rna_seq_reads, finite_ratios)
    #ax4.text(0.8, 0.9, "R %.2f" % pearson_r, transform=ax4.transAxes)
    #f.tight_layout()
    gs.update(wspace=.5, hspace=.5)
    print "\tWriting figure to %s..." % str(
        'figs/%s/%s_enrichment_ratio_scatterplot.pdf' % (top_level_dir, label))
    if not os.path.exists('figs/%s' % top_level_dir):
        os.system('mkdir ' + str('figs/%s' % top_level_dir))
    plt.savefig('figs/%s/%s_enrichment_ratio_scatterplot.pdf' % (
        top_level_dir, label))
    plt.clf()

def remove_all_but_highest_peak_per_gene(peaks):
    #gene_ranks = append(peaks.loc[:, ('gene_name', 'height')])
    peaks.sort('height', inplace=True, ascending=False)
    peaks.drop_duplicates('gene_name', inplace=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    Create a scatterplot comparing CLIP signal with a negativ IP and RNA-seq.
    Works on either FBF-1/FBF-2 separately, or together.''')
    parser.add_argument('-i', '--input', help='''Input directory.''')
    args = parser.parse_args()
    top_level_dir = args.input
    peaks = {}
#    rna_seq_bam_filename = '/home/dp/Desktop/bams/celegans/modencode_4594.bam'
    for filename in glob.glob(top_level_dir + '/combined*.txt'):
        print "\tLoading %s to make scatterplot..." % filename
        peaks[filename] = pandas.read_csv(filename, sep='\t')
        #peaks[filename] = peaks[filename].head(500)
        if re.search('fbf1', filename) is not None: label = 'fbf1'
        elif re.search('fbf2', filename) is not None: label = 'fbf2'
        else:
            print "Could not tell if filename %s was FBF-1 or FBF-2..." % filename
            sys.exit()
        peaks[filename] = rename_columns(peaks[filename], which=label)
        scatterplot_vs_controls(peaks[filename], label=label,
                                top_level_dir=top_level_dir)
        label += '_highest_peak_per_gene'
        remove_all_but_highest_peak_per_gene(peaks[filename])
        scatterplot_vs_controls(peaks[filename], label=label,
                                top_level_dir=top_level_dir)
