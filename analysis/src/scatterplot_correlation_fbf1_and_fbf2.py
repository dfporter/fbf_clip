import pandas
import re
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import matplotlib
import sys
from scipy import stats
import argparse

parser = argparse.ArgumentParser(
    description='''Scatterplot between fbf1 and fbf2.'''
)
parser.add_argument('-d', '--directory', dest='directory',
                    default='new_pypeaks_fdr1_negip_local/',
                    help='Input directory of peaks files.')

args = parser.parse_args()
top_level_dir = args.directory
peaks = []
gene_ranks = []
for filename in [top_level_dir + '/combined_fbf1.txt',
                 top_level_dir + '/combined_fbf2.txt']:
        peaks.append(pandas.DataFrame.from_csv(filename, sep='\t'))
        gene_ranks.append(peaks[-1].loc[:, ('gene_name', 'height')])
        gene_ranks[-1].sort('height', inplace=True, ascending=False)
        gene_ranks[-1].drop_duplicates('gene_name', inplace=True)
        gene_ranks[-1]['ranked'] = gene_ranks[-1]['height'].rank(ascending=False)
ranks = pandas.merge(gene_ranks[0], gene_ranks[1], on='gene_name', how='outer')
ranks['rank_dif'] = ranks['ranked_x'] - ranks['ranked_y']
ranks['height_dif'] = ranks['height_x'] - ranks['height_y']
fbf1 = peaks[0]
fbf2 = peaks[1]


matplotlib.rcParams.update({'font.size': 12})
plt.scatter(np.log10(ranks['height_x']), np.log10(ranks['height_y']),
            alpha=0.5, color='blue', s=4, edgecolor='none')
matplotlib.rcParams.update({'font.size': 6})
for index, row in ranks.iterrows():
    x = np.log10(ranks.loc[index, 'height_x'])
    y = np.log10(ranks.loc[index, 'height_y'])
    if y - x > np.log10(10):
        print "outlier: %s" % str((x,y))
        plt.annotate(ranks.loc[index, 'gene_name'], xy=(x,y),
                     textcoords='offset points', xytext = (-20, 20),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    if x - y > np.log10(10):
        print "outlier: %s" % str((x,y))
        plt.annotate(ranks.loc[index, 'gene_name'], xy=(x,y),
                     textcoords='offset points', xytext = (20, -20),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

plt.xlabel('log10 peak height FBF-1')
plt.ylabel('log10 peak height FBF-2')
outdir = 'figs/%s' % os.path.dirname(top_level_dir)
if not os.path.exists(outdir):
    os.system('mkdir ' + outdir)
plt.savefig(outdir + '/height_scatter_plot.pdf')
print "\tCalculating p values and rhos..."
(spearman_r, spearman_p_value) = stats.spearmanr(ranks.ranked_x, ranks.ranked_y)
ranks_in_both = pandas.merge(gene_ranks[0], gene_ranks[1], on='gene_name', how='inner')
(pearson_r, pearson_p_value) = stats.pearsonr(ranks_in_both.height_x, ranks_in_both.height_y)

table_dir = 'tables/%s' % os.path.basename(os.path.dirname(top_level_dir))
if not os.path.exists(table_dir):
    os.system('mkdir ' + table_dir)

with open(
    str(table_dir + '/fbf1_and_fbf2_correlation.txt'), 'w') as f:
    f.write("""
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
    pears_p_val=pearson_p_value))
