import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
#import pysam
import matplotlib
import sys
from matplotlib_venn import *
import os
import HTSeq
from subset_peaks_with_fbe import *
#from peak_locations import locate_in_gene, get_distance


def combine_peaks(replicates):
    combined = {}  # Key: (chrm, left, strand) tuple.
    first_rep_name = replicates.keys()[0]
    as_single_df = replicates[first_rep_name]
    as_single_df = pandas.concat([replicates[rep][:] for rep in replicates], ignore_index=True)
    print "after append: size of combined df: %i." % len(as_single_df)
    for index, row in as_single_df.iterrows():
        if index % 1000 == 0: print "combining %i/%i." % (index, len(as_single_df))
        gene = as_single_df.loc[index, 'gene_name']
        iv = [as_single_df.loc[index, 'chrm'], as_single_df.loc[index, 'left'],
              as_single_df.loc[index, 'right'], as_single_df.loc[index, 'strand']]
        overlapping_peaks = overlapping(iv, as_single_df)
        if len(overlapping_peaks) > 3:
            consensus = consensus_peak(overlapping_peaks)
            tup = (consensus['chrm'], consensus['left'], consensus['strand'])
            if tup not in combined:
                combined[tup] = consensus
        #if len(combined) > 10:
        #    return combined
    return combined


def combine_peaks_not_pandas(replicates):
    # Convert to a list of dicts.
    reps = {}
    for rep in replicates:
        reps[rep] = []
        for row in replicates[rep].itertuples():
            reps[rep].append(dict(
                (str(replicates[rep].columns[num-1]), val) for num, val in enumerate(row))
            )
    # Make a lookup by gene.
    by_gene = {}
    for rep in reps:
        by_gene[rep] = {}
        for row in reps[rep]:
            by_gene[rep].setdefault(row['gene_name'], [])
            by_gene[rep][row['gene_name']].append(row)
    # Get all genes.
    genes = []
    for rep in reps:
        genes.extend([x['gene_name'] for x in reps[rep]])
    genes = set(genes)
    overlaps = {}
    overlapping_peak_rows = {}
    for index, gene in enumerate(genes):
        _ranges = get_ranges_not_pandas(by_gene, gene)
        overlapping = get_overlapping_not_pandas(_ranges, gene)
        # Alters overlapping to remove shorter peaks that overlap
        # with higher peaks.
        remove_peaks_that_overlap_with_higher_peak(overlapping, _ranges)
        overlapping_peak_rows[gene] = []
        for peak_name in overlapping:
            overlapping_peak_rows[gene].append(_ranges[peak_name[0]][peak_name])
        #print c_peaks
    #print overlaps
    return overlapping_peak_rows


def consensus_peaks_not_pandas(overlapping, gene_ranges):
    # Do these peak ranges overlap?
    # Whenever two ranges overlap, take the highest.
    remove_peaks_that_overlap_with_higher_peak(overlapping, _ranges)
    return list_of_highest


def remove_peaks_that_overlap_with_higher_peak(overlapping, gene_ranges, verbose=False):
    still_overlapping = True
    to_remove = True
    while to_remove:
        to_remove = False
        for peak_a in overlapping:
            row_a = gene_ranges[peak_a[0]][peak_a]
            tup_a = (row_a['left'], row_a['right'])
            for peak_b in overlapping:
                row_b = gene_ranges[peak_b[0]][peak_b]
                tup_b = (row_b['left'], row_b['right'])
                if peak_a == peak_b: continue
                if overlap(tup_a, tup_b):
                    # Remove the shorter peak.
                    if row_a['height'] >= row_b['height']:
                        to_remove = peak_b
                    else:
                        to_remove = peak_a
                    still_overlapping = True
        if to_remove:
            #print "Removing %s" % str(to_remove)
            del overlapping[to_remove]
    if verbose:
        for peak_name in overlapping:
            print "%s " % str(peak_name)
            print "\t%i\t%i" % (gene_ranges[peak_name[0]][peak_name]['left'],
                                gene_ranges[peak_name[0]][peak_name]['right'])
            for overlapping_peak_name in overlapping[peak_name]:
                print "\tOverlaps %s" % str(overlapping_peak_name)
                print "\t%i\t%i" % (
                    gene_ranges[overlapping_peak_name[0]][overlapping_peak_name]['left'],
                    gene_ranges[overlapping_peak_name[0]][overlapping_peak_name]['right'])


def get_ranges_not_pandas(by_gene, gene):
    gene_ranges = {}
    for rep in by_gene:
        if gene not in by_gene[rep]:
            continue
        for peak in by_gene[rep][gene]:
            gene_ranges.setdefault(rep, {})
            peak_name = (rep, gene, peak['left'], peak['height'])
            gene_ranges[rep][peak_name] = peak
    return gene_ranges


def overlap(x, y):
    if x[0] <= y[0] <= x[1]:
        return True
    if x[0] <= y[1] <= x[1]:
        return True
    return False


def get_overlapping_not_pandas(_ranges, gene_name):
    """For each peak for this gene, find its overlapping peaks.
    Adds peak_name to overlaps for each peak for the gene,
    and assigns it to be a list of overlapping peaks (specifically,
    a list of their peak_names).
    Returns overlaps.
    """
    overlaps = {}
    for rep in _ranges:
        for peak_name in _ranges[rep]:
            row = _ranges[rep][peak_name]
            (num_rep_overlap, list_of_overlapping_peaks) = count_overlapping(
                row, _ranges)
            if num_rep_overlap >= 5:
                overlaps[peak_name] = list_of_overlapping_peaks
    return overlaps


def count_overlapping(cf_row, _ranges):
    """For a given peak row, returns the number of replicates
    with an overlapping peak and a list of peak_names of peaks
    overlapping the given peak.
    """
    list_of_overlapping_peaks = []
    cf_tup = (cf_row['left'], cf_row['right'])
    reps_with_overlapping_peak = set()
    for rep in _ranges:
        for peak_name in _ranges[rep]:
            row = _ranges[rep][peak_name]
            tup_b = (row['left'], row['right'])
            if overlap(cf_tup, tup_b):
                list_of_overlapping_peaks.append(peak_name)
                reps_with_overlapping_peak.add(rep)
    num_overlap = len(list(reps_with_overlapping_peak))
    return (num_overlap, list_of_overlapping_peaks)


def consensus_peak(peaks):
    rows = peaks[peaks['height']==max(peaks['height'])]
    return rows.iloc[0].to_dict()


def overlapping(iv, comb_df):
    overlap = comb_df[(comb_df['chrm']==iv[0]) & (comb_df['strand']==iv[3]) & (comb_df['left']<=iv[2]) & (iv[1]<=comb_df['right'])]
    return overlap


def write_combined_not_pandas(combined, filename, col_order):
    with open(filename, 'w') as f:
        header = ""
        for index, col in enumerate(col_order):
            if index > 0: header += "\t"
            header += str(col)
        f.write(header + "\n")
        all_rows = []
        for gene_name in combined:
            for peak_row in combined[gene_name]:
                all_rows.append(peak_row)
                #print peak_row
        all_rows = sorted(all_rows, key=lambda x: x['height'])[::-1]
        for peak_row in all_rows:
            li = ""
            for index, col in enumerate(col_order):
                if index > 0: li += "\t"
                if col in peak_row:
                    li += "%s" % str(peak_row[col])
            f.write(li + "\n")


if __name__ == '__main__':
    top_level_dir = sys.argv[1]
    replicates = {}
    col_order = []
    for filename in glob.glob(top_level_dir + '/*.txt'):
        print "\tLoading %s..." % filename
        replicates[filename] = pandas.read_csv(filename, sep='\t')
        col_order = replicates[filename].columns
    peak_rows_by_gene = combine_peaks_not_pandas(replicates)
    out_dir = top_level_dir.rstrip('/') + '_five_reps/'
    if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
    write_combined_not_pandas(peak_rows_by_gene, out_dir + 'combined_fbf.txt', col_order)
    sys.exit()
    old_version = '''
    combined_dict = combine_peaks(replicates)
    combined_df = pandas.DataFrame(combined.values())
    gtf_sep_cols = pandas.read_csv('lib/gtf_with_names_column.txt', sep='\t')
    locate_in_gene(gtf_sep_cols, combined_df)
    combined_df.to_csv('with_fbe_combined_fbfs_%s/%s' % (
        top_level_dir, 'fbf1_and_2'), sep='\t')
python src/filter_peaks_by_ratio_in_peak_region.py --analysis -i ${method} --output_dir filtered_${method}

    '''
