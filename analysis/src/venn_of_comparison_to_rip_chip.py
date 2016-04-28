import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib
import sys
from matplotlib_venn import *
import os

def make_venn_of_comparison_to_rip(peaks, label='no_label', subdir=''):
    targets = {}
    for exp in peaks:
        targets[exp] = set(peaks[exp]['gene_name'].dropna())
    make_venn(targets, label, subdir)

    
def make_venn(targets, label, subdir):
    plt.clf()
    overlaps = venn3([targets[exp] for exp in targets],
                     (exp for exp in targets))
    plt.savefig('figs/%s/overlap_with_rip_%s.pdf' % (
        subdir, label), format='pdf')
    plt.clf()
    overlap_fbf1 = venn2(
        [targets['RIP-ChIP'], targets['FBF-1']],
        ('RIP-ChIP', 'FBF-1'))
    plt.savefig('figs/%s/overlap_with_rip_%s.pdf' % (
        subdir, 'fbf1'), format='pdf')
    plt.clf()
    overlap_fbf2 = venn2(
        [targets['RIP-ChIP'], targets['FBF-2']],
        ('RIP-ChIP', 'FBF-2'))
    plt.savefig('figs/%s/overlap_with_rip_%s.pdf' % (
        subdir, 'fbf2'), format='pdf')
    plt.clf()

def run(top_dir):
    peaks = {}
    peaks['RIP-ChIP'] = pandas.read_csv('lib/ripchip/sd01.txt', sep='\t')
    peaks['RIP-ChIP'].columns = [
        'Affymetrix probe set ID', 'Gene', 'Ensembl Gene ID (unambigous only)',
        'WB Gene ID', 'gene_name', 'SAM score (d.value)', 'SAM rank',
        'Protein description', 'Protein remarks', """Annotated 3'UTR""",
        'No. of FBEs']
    peaks['RIP-ChIP'].drop_duplicates('gene_name', inplace=True)
    peaks['RIP-ChIP'] = peaks['RIP-ChIP'][0:1351]
    peaks['FBF-1'] = pandas.read_csv(top_dir + '/combined_fbf1.txt', sep='\t')
    peaks['FBF-2'] = pandas.read_csv(top_dir + '/combined_fbf2.txt', sep='\t')
    make_venn_of_comparison_to_rip(peaks,
                                   label='combined',
                                   subdir=os.path.dirname(top_dir))
    

if __name__ == '__main__':
    top_dir = sys.argv[1].rstrip('/') + '/'
    run(top_dir)    
