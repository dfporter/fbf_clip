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
import add_reads_columns

sys.path.insert(0, '/Network/Servers/file.biochem.wisc.edu/Volumes/BioSAN/Users/dfporter/.local/lib/python2.6/site-packages')
# Need HTSeq version 0.6, hence the sys.path alteration.
import HTSeq

#r = HTSeq.WiggleReader('fbf1/wigs/runx_fbf_gcca_20mapqper_mill_reads.wig')

# Apparently WiggleReader doesn't actually work for the bedgraph.
# So we write our own.

skip = r'''
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

'''
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
    arr2d = {ga.keys()[0]:[],
                          ga.keys()[1]:[]}
    for gene in peaks:
        reg = HTSeq.GenomicInterval(peaks[gene][0], peaks[gene][1], peaks[gene][2], peaks[gene][3])
        for exp in ga:
            integrals[exp].append(
                [gene, get_integral(ga[exp], reg)])
            arr2d[exp].append(integrals[exp][-1][1])
        if integrals[integrals.keys()[0]][-1][1] < 10 and integrals[integrals.keys()[1]][-1][1] < 10:
            print "{gene}:{iv}\nintegral {exp1}:{integ1}\n integral {exp2}:{integ2}\n".format(
                gene=gene, iv=str(reg),
                exp1=integrals.keys()[0], integ1=str(integrals[integrals.keys()[0]][-1]),
                exp2=integrals.keys()[1], integ2=str(integrals[integrals.keys()[1]][-1])
            )
    cutoff = 10
    for exp in integrals:
        integrals[exp] = filter(lambda x: x[1] > cutoff, integrals[exp])
    arr2d = [arr2d[arr2d.keys()[0]], arr2d[arr2d.keys()[1]]]
    print_stats(arr2d)
    return integrals

def print_stats(arr2d, top_level_dir=None):
    if top_level_dir is None: top_level_dir = 'fbf1_and_fbf2_correlation'
    (spearman_r, spearman_p_value) = stats.spearmanr(arr2d[0], arr2d[1])
    (pearson_r, pearson_p_value) = stats.pearsonr(arr2d[0], arr2d[1])

    if not os.path.exists('tables/%s' % os.path.dirname(top_level_dir)):
        os.system('mkdir ' + 'tables/%s' % os.path.dirname(top_level_dir))

    with open(
        str('tables/%s/' % os.path.dirname(top_level_dir)) + 'fbf1_and_fbf2_correlation.txt', 'w') as f:
        li = """
    Rank correlation, spearman rho: {s_rho}
    Rank correlation, spearman p-value for samples being the same: {spear_p_val}
    Peak height correlation for targets of both FBF-1 and FBF-2, Pearson R: {pears_r}
    \t...As R^2: {pears_r_squared}
    Two tailed p value for peak height correlation: {pears_p_val}
    """.format(
        s_rho=spearman_r,
        spear_p_val=spearman_p_value,
        pears_r=pearson_r,
        pears_r_squared=pearson_r**2,
        pears_p_val=pearson_p_value)
        print li
        f.write(li)


def get_integral(ga_obj, region):
    exp_score = 0.
    for iv, score in ga_obj[region].steps():
        exp_score += float(iv.end - iv.start) * score
    return exp_score

def convert_to_array(integrals):
    by_gene = convert_to_by_gene(integrals)
#    as_arr = []
    as_arr_with_gene = []
    for gene in by_gene:
        fbf1_val = [float(x[0]) for x in by_gene[gene] if x[1] == 'FBF-1'][0]
        fbf2_val = [float(x[0]) for x in by_gene[gene] if x[1] == 'FBF-2'][0]
#        as_arr.append(np.array([fbf1_val, fbf2_val]))
        as_arr_with_gene.append([fbf1_val, fbf2_val, gene])
    return as_arr_with_gene


def plot_integrals(as_arr, filename='integrals_in_peak_ranges.pdf',
                   input_dir='methods/', do_annotate=False):
    outliers_filename = 'data/%s/outliers_fbf1_vs_fbf2.txt' % os.path.basename(
        os.path.realpath(input_dir))
    r_row, r_p_value = sp.stats.pearsonr(
        [v[0] for v in as_arr], [v[1] for v in as_arr])  # Two-tailed p value
    spear_rho, spear_p_value = sp.stats.spearmanr(
        [v[0] for v in as_arr], [v[1] for v in as_arr])
    print """
Pearson R {pear}, p value {rpval}
Spearman R {spear_rho}, p value {spval}
    """.format(
        pear=r_row, rpval=r_p_value,
        spear_rho=spear_rho, spval=spear_p_value)
    plt.clf()
    plt.figure(figsize=(4,4))
    plt.axes().set_aspect('equal')
    plt.scatter([np.log10(v[0]) for v in as_arr],
                [np.log10(v[1]) for v in as_arr],
                alpha=0.5, color='black', s=4, edgecolor='none')
    max_val = np.max([[np.log10(v[0]) for v in as_arr],
                      [np.log10(v[1]) for v in as_arr]])
    matplotlib.rcParams.update({'font.size': 6})
    difs = {}
    for fbf1_val, fbf2_val, gene in as_arr:
        difs[gene] = (fbf1_val, fbf2_val,
                      np.log10(fbf1_val) - np.log10(fbf2_val), gene)
        if do_annotate:
            if fbf2_val/fbf1_val > 100.:
                if (fbf1_val > 1000) or (fbf2_val > 1000):
                    #print "%s %s" % (str(tup_x), str(tup_y))
                    plt.annotate(gene, xy=(np.log10(fbf1_val), np.log10(fbf2_val)),
                         textcoords='offset points', xytext = (-20, 20),
                         arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            if fbf1_val/fbf2_val > 100.:
                if (fbf1_val > 1000) or (fbf2_val > 1000):
                    plt.annotate(
                        gene, xy=(np.log10(fbf1_val), np.log10(fbf2_val)),
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
        f.write("\t".join(['FBF-1', 'FBF-2',
                'log10(FBF-1/FBF-2)', 'gene_name']) + '\n')
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
    plt.xlabel('$\mathregular{log_{10}}$ reads in peak FBF-1')
    plt.ylabel('$\mathregular{log_{10}}$ reads in peak FBF-2')
    plt.savefig(filename)
    plt.clf()


def run(args, ga):
    peaks = {}
    input_dir = args.directory
    for filename in glob.glob(input_dir + '/combined_fbf*.txt'):
        exp = os.path.basename(filename)
        peaks[exp] = pandas.read_csv(
            filename, sep='\t')
        peaks[exp] = peaks[exp][peaks[exp]['biotype']=='protein_coding']
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
    as_arr = convert_to_array(integrals) # list of (fbf1, fbf2, gene).
    get_two_fold_differences(as_arr, input_dir=input_dir)
    plot_integrals(as_arr, filename=plot_filename, input_dir=input_dir)

def convert_to_by_gene(integrals):
    exp1 = integrals.keys()[0]
    exp2 = integrals.keys()[1]
    if re.search('fbf2', exp1) is not None:
        x_label = 'FBF-2'
        y_label = 'FBF-1'
    else:
        x_label = 'FBF-1'
        y_label = 'FBF-2'
    by_gene = {}
    for tupx in integrals[exp1]:
        by_gene.setdefault(tupx[0], [])
        by_gene[tupx[0]].append((tupx[1], x_label))
    for tupx in integrals[exp2]:
        by_gene.setdefault(tupx[0], [])
        by_gene[tupx[0]].append((tupx[1], y_label))
    by_gene = only_joint_targs(by_gene)
    return by_gene


def get_two_fold_differences(as_arr, input_dir='methods.'):
    twofold_fbf1_up = 'data/%s/fbf1_two_fold_up.txt' % os.path.basename(
        os.path.realpath(input_dir))
    twofold_fbf2_up = 'data/%s/fbf2_two_fold_up.txt' % os.path.basename(
        os.path.realpath(input_dir))
    ratios = {}
    for fbf1, fbf2, gene in as_arr:
        ratios[gene] = [float(fbf1)/float(max([1., fbf2])), fbf1, fbf2]
    fbf1_up = dict((key, value) for key, value in ratios.items() if value[0]>2.)
    fbf2_up = dict((key, value) for key, value in ratios.items() if value[0]<0.5)


    fbf1_up_keys = sorted(fbf1_up.keys(), key=lambda x: fbf1_up[x][0], reverse=True)
    fbf2_up_keys = sorted(fbf2_up.keys(), key=lambda x: fbf2_up[x][0], reverse=False)
    with open(twofold_fbf1_up, 'w') as f:
        f.write('gene_name\tfbf1_reads\tfbf2_reads\rratio (FBF-1/FBF-2)\n')
        f.write("\n".join(
            ["{g}\t{f1}\t{f2}\t{r}".format(
                g=gene, f1=ratios[gene][1], f2=ratios[gene][2], r=ratios[gene][0]
            ) for gene in fbf1_up_keys]
        ))
    with open(twofold_fbf2_up, 'w') as f:
        f.write('gene_name\tfbf1_reads\tfbf2_reads\tratio (FBF-2/FBF-1)\n')
        f.write("\n".join(
            ["{g}\t{f1}\t{f2}\t{r}".format(
                g=gene, f1=ratios[gene][1], f2=ratios[gene][2], r=ratios[gene][0]**-1
            ) for gene in fbf2_up_keys]
        ))

def only_joint_targs(by_gene):
    f_by_gene = {}
    for gene in by_gene:
        if not (('FBF-1' in [x[1] for x in by_gene[gene]]) and (
            'FBF-2' in [x[1] for x in by_gene[gene]]
        )):
            continue
        f_by_gene[gene] = by_gene[gene]
    return f_by_gene


def get_ga_from_bedgraph(args, ga=None, bed_type='coverage'):
    if ga is None: ga = {}
    if args.load_data_from_pickle:
        with open('data/integrals_from_scatterplot_correlation_by_wig.p', 'r') as f:
            integrals = pickle.load(f)
    else:
        ga['combined_fbf1.txt'] = add_reads_columns.load_bedgraph(
            'bedgraph_norm/combined_fbf1.wig')
        ga['combined_fbf2.txt'] = add_reads_columns.load_bedgraph(
            'bedgraph_norm/combined_fbf2.wig')
        skip = r'''
        for filename_list in [
            ('bedgraph_norm/combined_fbf1_+.wig',
             'bedgraph_norm/combined_fbf1_-.wig',
             'combined_fbf1.txt'),
            ('data/wigs_five_prime/fbf2_reads_plus.bed',
             'data/wigs_five_prime/fbf2_reads_minus.bed',
             'combined_fbf2.txt')]:
            peaks_filename = filename_list[2]
            # Build the ga dict with keys by peaks_fileame
            #load_bedgraph(filename_list, ga, use_key=peaks_filename)'''
    return ga


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
#    parser.add_argument('-s', '--cols', default=False, action='store_true',
#                        help='''Use the reads_fbf1/reads_fbf2 columns from the \
#                        input peaks file, as generated by \
#                        filter_peaks_by_ratio_in_peak_region.py --load_bed_info''')
    args = parser.parse_args()
    ga = {}
    get_ga_from_bedgraph(args, ga=ga)
    run(args, ga)

