import sys
import glob
import os
import re
import pandas
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
import pickle
import argparse

sys.path.insert(0, '/Network/Servers/file.biochem.wisc.edu/Volumes/BioSAN/Users/dfporter/.local/lib/python2.6/site-packages')
# Need HTSeq version 0.6, hence the sys.path alteration.
import HTSeq

#r = HTSeq.WiggleReader('fbf1/wigs/runx_fbf_gcca_20mapqper_mill_reads.wig')

# Apparently WiggleReader doesn't actually work for the bedgraph.
# So we write our own.

def load_bedgraph(filename_list, ga, use_key=False):
    if use_key:
        exp = use_key
    else: exp = filename_list[0]
    ga[exp] = HTSeq.GenomicArray(chroms='auto', stranded=True)
    ratio_fbf1_to_2 = float(9792191+3166675+10408265)/float(7680463+884888+5584323)
    #rep = HTSeq.WiggleReader(filename)
    if exp == 'combined_fbf2.txt':
        with open(filename_list[0], 'r') as f:
            next(f)
            for line in f:
                s = line.rstrip('\n').split('\t')
                ga[exp][HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] = ratio_fbf1_to_2 * float(s[3])
        with open(filename_list[1], 'r') as f:
            next(f)
            for line in f:
                s = line.rstrip('\n').split('\t')
                ga[exp][HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] = ratio_fbf1_to_2 * float(s[3])
    else:
        with open(filename_list[0], 'r') as f:
            next(f)
            for line in f:
                s = line.rstrip('\n').split('\t')
                ga[exp][HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] = float(s[3])
        with open(filename_list[1], 'r') as f:
            next(f)
            for line in f:
                s = line.rstrip('\n').split('\t')
                ga[exp][HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] = float(s[3])


def get_all_ranges(peaks):
    by_rep_by_gene = {}
    for filename in peaks:
        by_rep_by_gene[filename] = {}
        for index, row in peaks[filename].iterrows():
            iv = (row['chrm'], int(row['left']), int(row['right']),
                  row['strand'], float(row['height']))
            by_rep_by_gene[filename].setdefault(row['gene_name'], [])
            by_rep_by_gene[filename][row['gene_name']].append(iv)
    exp1 = peaks.keys()[0]
    exp2 = peaks.keys()[1]
    genes1 = set(by_rep_by_gene[exp1].keys())
    genes2 = set(by_rep_by_gene[exp2].keys())
    highest = {}
    all_genes = genes1 | genes2
    for gene in all_genes:
        if gene in by_rep_by_gene[exp1]:
            peaks1 = by_rep_by_gene[exp1][gene]
        else: peaks1 = []
        if gene in by_rep_by_gene[exp2]:
            peaks2 = by_rep_by_gene[exp2][gene]
        else: peaks2 = []
        highest[gene] = highest_peak_in_gene(peaks1, peaks2)
    return highest


def highest_peak_in_gene(peaks1, peaks2):
    peaks1.extend(peaks2)
    return sorted(peaks1, key=lambda x: float(x[4]))[-1]


def highest_overlapping_one_peak_from_a_list(peak1, peaks2):
    """Return a list of peak
    """
    overlapping_with_peak1 = set([peak1])
    for peak2 in peaks2:
        overlap = False
        if peak1[1] <= peak2[1] <= peak1[2]:
            overlap = True
        if peak1[1] <= peak2[2] <= peak1[2]:
            overlap = True
        if overlap:
            overlapping_with_peak1.add(peak2)
    overlapping_with_peak1 = list(overlapping_with_peak1)
    highest_overlapping.add(
        sorted(
            overlapping_with_peak1, key=lambda x: x[4])[-1])
    return highest_overlapping


def get_signal_in_ranges(peaks, ga):
    integrals = {ga.keys()[0]: [],
                 ga.keys()[1]: []}
    for gene in peaks:
        reg = HTSeq.GenomicInterval(peaks[gene][0], peaks[gene][1], peaks[gene][2], peaks[gene][3])
        if gene == 'pyr-1':
            print "get_signal_in_ranges, pyr-1: reg={reg}".format(reg=str(reg))
            for exp in ga:
                print str(exp)
                total = 0.
                for iv, score in ga[exp][reg].steps():
                    print iv, score
                    total += float(iv.end - iv.start) * score
                    print "total so far=%f" % total
        for exp in ga:
            integrals[exp].append(
                [gene, get_integral(ga[exp], reg)])
        if integrals[integrals.keys()[0]][-1][1] < 10 and integrals[integrals.keys()[1]][-1][1] < 10:
            print "{gene}:{iv}\nintegral {exp1}:{integ1}\n integral {exp2}:{integ2}\n".format(
                gene=gene, iv=str(reg),
                exp1=integrals.keys()[0], integ1=str(integrals[integrals.keys()[0]][-1]),
                exp2=integrals.keys()[1], integ2=str(integrals[integrals.keys()[1]][-1])
            )
    return integrals


def get_integral(ga_obj, region):
    exp_score = 0.
    for iv, score in ga_obj[region].steps():
        exp_score += float(iv.end - iv.start) * score
    return exp_score


def plot_integrals(integrals, filename='integrals_in_peak_ranges.pdf',
                   outliers_filename='outliers_fbf1_vs_fbf2.txt'):
    exp1 = integrals.keys()[0]
    exp2 = integrals.keys()[1]
    if re.search('fbf2', exp1) is not None:
        x_label = 'FBF-2'
        y_label = 'FBF-1'
    else:
        x_label = 'FBF-1'
        y_label = 'FBF-2'
    x = sp.array([tup[1] for tup in integrals[exp1]])
    y = sp.array([tup[1] for tup in integrals[exp2]])
    r_row, r_p_value = sp.stats.pearsonr(x, y)  # Two-tailed p value
    spear_rho, spear_p_value = sp.stats.spearmanr(x, y)
    print """
Pearson R {pear}, p value {rpval}
Spearman R {spear_rho}, p value {spval}
    """.format(
        pear=r_row, rpval=r_p_value,
        spear_rho=spear_rho, spval=spear_p_value)
    plt.clf()
    plt.figure(figsize=(4,4))
    plt.axes().set_aspect('equal')
    plt.scatter(np.log10(x), np.log10(y),
                alpha=0.5, color='black', s=4, edgecolor='none')
    max_val = np.max([np.log10(x), np.log10(y)])
    matplotlib.rcParams.update({'font.size': 6})
    difs = {}
    for index, tup_x in enumerate(integrals[exp1]):
        gene = tup_x[0]
        tup_y = integrals[exp2][index]
        if tup_x[0] != tup_y[0]:
            print "Genes not equal? %s %s" % (str(tup_x), str(tup_y))
        difs[gene] = (tup_x[1], tup_y[1], np.log10(tup_y[1]) - np.log10(tup_x[1]), gene)
        if np.log10(tup_y[1]) - np.log10(tup_x[1]) > np.log10(4):
            if (tup_x[1] > 1000) or (tup_y[1] > 1000):
                #print "%s %s" % (str(tup_x), str(tup_y))
                plt.annotate(gene, xy=(np.log10(tup_x[1]), np.log10(tup_y[1])),
                     textcoords='offset points', xytext = (-20, 20),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        if np.log10(tup_x[1]) - np.log10(tup_y[1]) > np.log10(16):
            if (tup_x[1] > 1000) or (tup_y[1] > 1000):
                #print "%s %s" % (str(tup_x), str(tup_y))
                plt.annotate(gene, xy=(np.log10(tup_x[1]), np.log10(tup_y[1])),
                     textcoords='offset points', xytext = (20, -20),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    # Add lines at a given fold change difference.
    fold_cutoff = 2
    offset = np.log10([fold_cutoff])[0]
    # plt.plot([x1, x2], [y1, y2]) to plot (x1, y1) to (x2, y2).
    plt.plot([0, 1e6], [0 + offset, 1e6 + offset], 'r--')
    plt.plot([0, 1e6], [0 - offset, 1e6 - offset], 'r--')
    difs_list = sorted(difs.values(), key=lambda x: x[2])[::-1]
    f_1_up = open('data/fbf1_enriched_targets.txt', 'w')
    f_2_up = open('data/fbf2_enriched_targets.txt', 'w')
    with open(outliers_filename, 'w') as f:
        f.write("\t".join([x_label, y_label,
                'log10(' + y_label + '/' + x_label + ')', 'gene_name']) + '\n')
        for tup in difs_list:
            if (tup[0] > 600) or (tup[1] > 600):
                if np.log10(tup[1]) - np.log10(tup[0]) > np.log10(2):
                    f_1_up.write("%s\n" % tup[-1])
                if np.log10(tup[0]) - np.log10(tup[1]) > np.log10(2):
                    f_2_up.write("%s\n" % tup[-1])
                f.write("\t".join([str(x) for x in tup]) + '\n')
    f_1_up.close()
    f_2_up.close()
    plt.xlim([0, np.round(max_val, decimals=0)])
    plt.ylim([0, np.round(max_val, decimals=0)])
    plt.xlabel('$\mathregular{log_{10}}$ reads in peak ' + x_label)
    plt.ylabel('$\mathregular{log_{10}}$ reads in peak ' + y_label)
    plt.savefig(filename)
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Add sequence, has_fbe, rip_column, ect.
        info columns to a given directory -d of peaks files.
        Can also combine replicates and output the combined list.
        Will apply a normalization ratio of ~1.6 to FBF-2 reads.'''
    )
    parser.add_argument('-d', '--directory', dest='directory',
                        default='new_pypeaks_fdr1_negip_local/',
                        help='Input directory of peaks files.')
    parser.add_argument('-p', '--load_data_from_pickle',
                        default=False,
                        action='store_true',
                        help='''Load data from the file \
                        data/integrals_from_scatterplot_correlation_by_wig.p''')
#    parser.add_argument('-c', '--cols', default=False, action='store_true',
#                        help='''Use the reads_fbf1/reads_fbf2 columns from the \
#                        input peaks file, as generated by \
#                        filter_peaks_by_ratio_in_peak_region.py --load_bed_info''')
    args = parser.parse_args()
    input_dir = args.directory
    #globstr = input_dir + '/*per_mill_reads.wig'
    ga = {}
    if args.load_data_from_pickle:
        with open('data/integrals_from_scatterplot_correlation_by_wig.p', 'r') as f:
            integrals = pickle.load(f)
    else:
        for filename_list in [
            ('data/wigs_five_prime/fbf1_reads_plus.bed',
             'data/wigs_five_prime/fbf1_reads_minus.bed',
             'combined_fbf1.txt'),
            ('data/wigs_five_prime/fbf2_reads_plus.bed',
             'data/wigs_five_prime/fbf2_reads_minus.bed',
             'combined_fbf2.txt')]:
            peaks_filename = filename_list[2]
            # Build the ga dict with keys by peaks_fileame
            load_bedgraph(filename_list, ga, use_key=peaks_filename)
        peaks = {}
        for filename in glob.glob(input_dir + '/combined_fbf*.txt'):
            peaks[os.path.basename(filename)] = pandas.DataFrame.from_csv(filename, sep='\t')
        highest_per_gene = get_all_ranges(peaks)
        integrals = get_signal_in_ranges(highest_per_gene, ga)
        with open('data/integrals_from_scatterplot_correlation_by_wig.p', 'w') as f:
            pickle.dump(integrals, f)
    plot_filename = 'figs/%s/scatterplot_by_wig.pdf' % os.path.basename(
        os.path.realpath(input_dir))
    if not os.path.exists('data/%s/' % os.path.basename(
        os.path.realpath(input_dir))):
        os.system('mkdir ' + 'data/%s/' % os.path.basename(
        os.path.realpath(input_dir)))
    outliers_filename = 'data/%s/outliers_fbf1_vs_fbf2.txt' % os.path.basename(
        os.path.realpath(input_dir))
    plot_integrals(integrals, filename=plot_filename,
                   outliers_filename=outliers_filename)
