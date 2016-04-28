import pandas
import csv
import re


def get_ripchip_targets(get_sam_rank=False):
    ripchip_filename = 'lib/ripchip/sd01.txt'
    ripchip = pandas.read_csv(
        ripchip_filename, sep='\t')
    top_1350 = []
    by_sam_rank = {}
    for index, row in ripchip.iterrows():
        if row['Gene public name'] not in top_1350:
            top_1350.append(row['Gene public name'])
            by_sam_rank[row['Gene public name']] = row['SAM rank']
        if len(list(set(top_1350))) >= 1350:
            break
    #top_1350_rows = ripchip[0:1349]  # There is a row zero.
    #rip_targets = list(top_1350['Gene public name'])
    print "Top targets: {t}.\nLast targets: {b}.".format(
        t=str(top_1350[0:2]), b=str(top_1350[-2:]))
    if get_sam_rank:
        return top_1350, by_sam_rank
    return top_1350


def add_column_of_overlap(df, rip_targets):
    rip_targets = set(rip_targets)
    df['is_rip_chip_target'] = 0
    for index, row in df.iterrows():
        if df.loc[index, 'gene_name'] in rip_targets:
            df.loc[index, 'is_rip_chip_target'] = 1
    
if __name__ == '__main__':
    import sys
    try:
        fname = sys.argv[1]
    except:
        print "Pass a peaks filename as argument."
        sys.exit()
    import pandas
    peaks = pandas.read_csv(fname, sep='\t', index_col=None)
    rip, sam_ranks = get_ripchip_targets(get_sam_rank=True)
    num_clip = len(list(set(peaks['gene_name'].tolist())))
    num_rip = len(list(set(rip)))
    num_only_clip =len(set(peaks.gene_name) - set(rip))
    num_only_rip = len(set(rip) - set(peaks.gene_name))
    num_both = len(
        set(rip) & set(peaks.gene_name))
    num_either = len(set(rip) | set(peaks.gene_name))
    num_neither = 20000 - num_either
    fraction_of_clip = float(num_both)/float(num_clip)
    fraction_of_rip = float(num_both)/float(num_rip)
    table = [
        [num_only_clip, num_both],
        [num_neither, num_only_rip]]
    print '''
       [[num_only_clip, num_both],
        [num_neither, num_only_rip]]'''
    print table
    print "New targets: %i" % num_only_clip
    print "Fraction of CLIP in RIP: %f" % float(fraction_of_clip)
    print "Fraction of RIP in CLIP: %f" % float(fraction_of_rip)

    import scipy
    import scipy.stats
    res = scipy.stats.fisher_exact(table)
    print res
    max_heights = {}
    peaks_d = peaks.to_dict('records')
    peak_heights = {}
    for row in peaks_d:
        peak_heights[row['gene_name']] = row['height']
    arr2d = [[], []]
    for gene in list(set(rip) & set(peaks['gene_name'].tolist())):
        arr2d[0].append(sam_ranks[gene])
        arr2d[1].append(peak_heights[gene])
    print scipy.stats.spearmanr(arr2d[0], arr2d[1])
