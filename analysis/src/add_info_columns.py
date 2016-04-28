"""

Adds several columns of information to input files.
Outputs to a new directory or edits the input file.
Unclear which is better.

"""
import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
#from scipy.optimize import minimize_scalar
#import pysam
import matplotlib
import sys
#from matplotlib_venn import *
import os
import HTSeq
import compare_with_ripchip
import subset_peaks_with_fbe
import argparse
#from peak_locations import locate_in_gene, get_distance
import time

def df_to_dict(df, key='gene_name'):
    d = {}
    for index, row in df.iterrows():
        d.setdefault(row[key], [])
        d[row[key]].append(row.to_dict())
    return d

# Duplicated from peak_locations.py.
def locate_in_gene(gtf_sep_cols, combined, use_this_column=None):
    print "\tLocating peaks..."
    start_time = time.time()
    if type(gtf_sep_cols) == pandas.DataFrame:
        gtf_sep_cols = df_to_dict(gtf_sep_cols)
        print "\tConverted to dict: time %f m" % (float(time.time() - start_time)/60.)
        start_time = time.time()
    combined['biotype'] = 'Unknown'
    if use_this_column is not None:
        if use_this_column in combined.columns:
            print "Using %s to set the peak location." % use_this_column
        else:
            print "Asked to use the column %s to set peak location, but \
was not found... Reverting to use where the majority of the \
peak region lies." % use_this_column
            use_this_column = None
    for index, peak_row in combined.iterrows():
        gene = str(combined.loc[index, 'gene_name'])
#        rows = gtf_sep_cols[gtf_sep_cols['gene_name']==gene]
        rows = gtf_sep_cols[gene]
        if len(rows) == 0:
            combined.loc[index, 'location'] = "Unknown"
            continue
        combined.loc[index, 'biotype'] = rows[0]['biotype']
        rows = [x for x in rows if x['2']=='CDS']
        if len(rows) == 0:
            combined.loc[index, 'location'] = "ncRNA"
            continue
        gene_left = min([x['3'] for x in rows])
        gene_right = max([x['4'] for x in rows])
        #gene_strand = rows['6']
        if use_this_column is not None:
            pos = peak_row[use_this_column]
            if (pos - gene_right) > 0:
                if combined.loc[index, 'strand'] == '+':
                    combined.loc[index, 'location'] = "3'UTR"
                if combined.loc[index, 'strand'] == '-':
                    combined.loc[index, 'location'] = "5'UTR"
            elif (pos - gene_left) < 0:
                if combined.loc[index, 'strand'] == '+':
                    combined.loc[index, 'location'] = "5'UTR"
                if combined.loc[index, 'strand'] == '-':
                    combined.loc[index, 'location'] = "3'UTR"
            else:
                combined.loc[index, 'location'] = "CDS"
            continue
        dist = get_distance((gene_left, gene_right),
                            (combined.loc[index, 'left'],
                             combined.loc[index, 'right']))
        if combined.loc[index, 'strand'] == '-':
            dist = -1 * dist
        combined.loc[index, 'dist_to_CDS'] = dist
        if not dist:
            combined.loc[index, 'location'] = "CDS"
        if dist < 0:
            combined.loc[index, 'location'] = '''5'UTR'''
        if dist > 0:
            combined.loc[index, 'location'] = '''3'UTR'''
    print "\tEdited peaks: time %f m" % (float(time.time() - start_time)/60.)


# Duplicated from peak_locations.py
def get_distance(geneiv, peakiv):
    if geneiv[1] < peakiv[0]:
        return peakiv[0] - geneiv[1]
    if peakiv[1] < geneiv[0]:
        return peakiv[1] - geneiv[0]
    if geneiv[0] < peakiv[0] < geneiv[1]:
        if peakiv[1] < geneiv[1]:
            return 0  # Entirely in CDS.
        else:
            in_cds = geneiv[1] - peakiv[0]
            frac_in_cds = float(in_cds)/float(peakiv[1] - peakiv[0])
            if frac_in_cds > 0.5:
                return 0
            else:
                midpoint = float(peakiv[1] + peakiv[0])/2.0
                return  midpoint - geneiv[1]
    if geneiv[0] < peakiv[1] < geneiv[1]:
        if peakiv[0] > geneiv[0]:
            return 0  # Entirely in CDS.
        else:
            in_cds = peakiv[1] - geneiv[0]
            frac_in_cds = float(in_cds)/float(peakiv[1] - peakiv[0])
            if frac_in_cds > 0.5:
                return 0
            else:
                midpoint = float(peakiv[1] + peakiv[0])/2.0
                return  midpoint - geneiv[0]
    return 0


def combine_peaks(replicates):
#    print 'passed %s' % str(replicates)
    combined = {}  # Key: (chrm, left, strand) tuple.
    first_rep_name = replicates.keys()[0]
    as_single_df = replicates[first_rep_name]
    as_single_df = pandas.concat([replicates[rep][:] for rep in replicates], ignore_index=True)
    print "After append: size of combined df: %i." % len(as_single_df)
    for index, row in as_single_df.iterrows():
        if index % 2000 == 1:
            print "Combining %i/%i." % (index, len(as_single_df))
        gene = as_single_df.loc[index, 'gene_name']
        iv = [as_single_df.loc[index, 'chrm'], as_single_df.loc[index, 'left'],
              as_single_df.loc[index, 'right'], as_single_df.loc[index, 'strand']]
        overlapping_peaks = overlapping(iv, as_single_df)
        if len(overlapping_peaks) > 4:
            consensus_dict, consensus_row = consensus_peak(overlapping_peaks)
            tup = (consensus_dict['chrm'],
                   consensus_dict['left'],
                   consensus_dict['strand'])
            if tup not in combined:
                combined[tup] = consensus_dict
        #if len(combined) > 10:
        #    return combined
    combined_df = pandas.DataFrame.from_dict(combined, orient='index')
    print combined_df.head()
    return combined_df


def consensus_peak(peaks):
    rows = peaks[peaks['height']==max(peaks['height'])]
    return (rows.iloc[0].to_dict(), rows.iloc[0])


def overlapping(iv, comb_df):
    overlap = comb_df[(comb_df['chrm']==iv[0]) & (comb_df['strand']==iv[3]) & (comb_df['left']<=iv[2]) & (comb_df['right']>=iv[1])]
    return overlap


def add_info(peaks, rip_targets, gtf_sep_cols,
             gtf_sep_cols_dict, peaks_fname, sequences,
             deseq=None, cufflinks=None,
             ):
    """Adds sequences, fbes, RIP-chip and peak location.
    Takes a dataframe.
    Returns the altered dataframe."""
    if (deseq is not None) and (cufflinks is not None):
    # Compare with the list in Ortiz et al. of gonadal genes.
        add_gonad_expression(peaks, deseq, cufflinks)
    # Add genomic sequence in the peak.
    peaks = subset_peaks_with_fbe.add_seqs(peaks, sequences)
    # Get yes/no FBE, -1/-2 C, number of FBEs.
    subset_peaks_with_fbe.score_binding_site(peaks)
    # Overlaps with the top 1350 unique RIP-chip targets?
    compare_with_ripchip.add_column_of_overlap(
        peaks, rip_targets)
    # Locate each peak in the gene.
    if re.search('fbf1', peaks_fname):
        locate_in_gene(
            gtf_sep_cols_dict, peaks, use_this_column='fbf1_reads_pos_of_max_coverage')
    elif re.search('fbf2', peaks_fname):
        locate_in_gene(
            gtf_sep_cols_dict, peaks, use_this_column='fbf2_reads_pos_of_max_coverage')
    else:
        locate_in_gene(
            gtf_sep_cols_dict, peaks)
    return peaks


def add_gonad_expression(peaks, ds, cf):
    for index, row in peaks.iterrows():
        if peaks.loc[index, 'gene_name'] in ds:
            peaks.loc[index, 'Present_in_gonad_by_DESeq'] = 1
        else:
            peaks.loc[index, 'Present_in_gonad_by_DESeq'] = 0
        if peaks.loc[index, 'gene_name'] in cf:
            peaks.loc[index, 'Present_in_gonad_by_Cufflinks'] = 1
        else:
            peaks.loc[index, 'Present_in_gonad_by_Cufflinks'] = 0
    return peaks


def read_deseq_file(ds_fname='lib/ortiz/DESeq_genes_in_gonad.txt'):
    ds = {}
    with open(ds_fname, 'r') as f:
        header = next(f)
        header = header.rstrip('\n').split('\t')
        for li in f:
            s = li.rstrip('\n').split('\t')
            ds[s[0]] = dict((header[i], s[i]) for i in range(0, len(s)))
    return ds

def read_cufflinks_file(ds_fname='lib/ortiz/Cufflinks_genes_in_gonad.txt'):
    ds = {}
    with open(ds_fname, 'r') as f:
        header = next(f)
        header = header.rstrip('\n').split('\t')
        for li in f:
            s = li.rstrip('\n').split('\t')
            ds[s[0]] = dict((header[i], s[i]) for i in range(0, len(s)))
    return ds


def get_columns(df):
    obs_col = set([str(x) for x in df.columns])
    exp_col = set(minimal_columns)
    found_col = obs_col & exp_col
    missing_col = exp_col - obs_col
    found_cols_to_output = [x for x in minimal_columns if x in found_col]
    print "****\nget_columns:"
    print found_col
    print missing_col
    print '****'
    return found_cols_to_output


def add_info_to_directory(
        args, rip_targets, gtf_sep_cols, deseq, cufflinks, minimal_columns):
    # Add info columns to the peaks in the input dir and either
    # edit in place or output to the given directory.
    peaks_d = {}
    start_time = time.time()
    if not args.skip_location:
        gtf_sep_cols_dict = df_to_dict(gtf_sep_cols)
    fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    print "\tConverted to dict: time %f m" % (float(time.time() - start_time)/60.)
    for peaks_fname in glob.glob(args.directory + '/*.txt'):
        print ">%s" % peaks_fname
        if len(open(peaks_fname).readlines()) < 2: continue
        peaks_d[peaks_fname] = pandas.read_csv(peaks_fname, sep='\t')
        print "...Loaded file"
        if args.skip_location:
            peaks_d[peaks_fname] = add_minimal_info(peaks_d[peaks_fname],
                    rip_targets, peaks_fname, sequences)
        else:
            peaks_d[peaks_fname] = add_info(peaks_d[peaks_fname], rip_targets,
                                        gtf_sep_cols, gtf_sep_cols_dict,
                                        peaks_fname, deseq, cufflinks,
                                        sequences)
    for filename in [x for x in peaks_d if peaks_d[x] is not None]:
        if args.fewer_columns:
            found_cols_to_output = get_columns(peaks_d[filename])
            peaks_d[filename] = peaks_d[filename][found_cols_to_output]
        if args.edit_in_place: outfname = filename
        else:
            outfname = '%s/%s' % (args.output_dir, os.path.basename(filename))
        if args.fewer_columns:
            peaks_d[filename].to_csv(outfname, sep='\t', index=False,
                                     columns=found_cols_to_output)
        else:
            peaks_d[filename].to_csv(outfname,  sep='\t', index=False)


def add_minimal_info(peaks, rip_targets,
                    peaks_fname, sequences):
        # Compare with the list in Ortiz et al. of gonadal genes.
        # Add genomic sequence in the peak.
        print '1'
        peaks = subset_peaks_with_fbe.add_seqs(peaks, sequences)
        # Get yes/no FBE, -1/-2 C, number of FBEs.
        print '2'
        subset_peaks_with_fbe.score_binding_site(peaks)
        # Overlaps with the top 1350 unique RIP-chip targets?
        print '3'
        compare_with_ripchip.add_column_of_overlap(
            peaks, rip_targets)
        print 'done'

if __name__ == '__main__':
    #print '\tTurn on locate_in_gene again!'
    parser = argparse.ArgumentParser(
        description='''Add sequence, has_fbe, rip_column, ect. \
        info columns to a given directory -d of peaks files. \
        Can also combine replicates and output the combined list.'''
    )
    parser.add_argument('-d', '--directory', dest='directory',
                        default='new_pypeaks_fdr1_negip_local/',
                        help='Input directory of peaks files.')
    parser.add_argument('-e', '--edit_in_place', default=False,
                        action='store_true', help='''Add info columns \
                        to the peaks files in-place. \
                        Default: Output to a with_info_* directory.''')
    parser.add_argument('-c', '--combine_replicates',
                        default=False, action='store_true',
                        help='''Combine peaks in 5/6 replicates.''')
    parser.add_argument('-s', '--skip_location', action='store_true',
                        default=False)
    parser.add_argument('-o', '--output_dir',
                        default=None, help='''Output directory for a \
                        peaks with info and a combined list of 5/6 replicate peaks. \
                        Default: with_info_<input_dir>/''')
    parser.add_argument('-w', '--fewer_columns', default=False,
                        action='store_true', help='''Reduce the number \
                        of info columns in the output.''')
    args = parser.parse_args()
    args.directory = args.directory.rstrip('/') + '/'
    # Setup the output directory in case we use it.
    if args.edit_in_place:
        args.output_dir = args.directory
    if (args.output_dir is None):
        args.output_dir = 'with_info_%s' % os.path.basename(os.path.dirname(args.directory))
    if not os.path.exists(args.output_dir):
        os.system('mkdir ' + args.output_dir)
    # Get some info from library files.
    rip_targets = compare_with_ripchip.get_ripchip_targets()
    gtf_sep_cols = pandas.read_csv('lib/gtf_with_names_column.txt', sep='\t')
    deseq = read_deseq_file()
    cufflinks = read_cufflinks_file()
    # Read in all the peaks in the given directory.
    # If asked with --combine_replicates, combine replicates and write to --output_combined_dir.
    global minimal_columns
    minimal_columns = [
            '0', '1', '2', '3', '4', '5', '6', '7', '8',
            'chrm', 'gene_name', 'height', 'left', 'right',
            'strand', 'transcript_id', 'transcript_name',
            'seq', 'has_fbe', 'number_of_fbes',
            'minus_one_c', 'minus_two_c', 'minus_one_or_two_c',
            'is_rip_chip_target', 'dist_to_CDS', 'location',
            'neg_ip_local_nb', 'neg_ip_local_norm', 'neg_ip_gene_nb',
            'neg_ip_gene_norm',
            'fbf1_reads', 'fbf2_reads', 'n2_from_fbf1_reads',
            'n2_from_fbf2_reads', 'rna_seq_reads',
            'biotype',]
    add_info_to_directory(
        args, rip_targets, gtf_sep_cols, deseq, cufflinks, minimal_columns)
    if args.combine_replicates:
        replicates = {}
        for filename in glob.glob(args.directory + '/fbf*'):
            replicates[filename] = pandas.read_csv(filename, sep='\t')
        combined_df = combine_peaks(replicates)
        five_reps_dir = args.directory
        five_reps_dir = five_reps_dir.rstrip('/') + '_five_reps/'
        os.system('mkdir ' + five_reps_dir)
        if args.fewer_columns:
            reduced_col_df = combined_df[get_columns(combined_df)]
            reduced_col_df.to_csv('%s/combined_fbf.txt' % five_reps_dir, sep='\t',)
        else:
            combined_df.to_csv('%s/combined_fbf.txt' % five_reps_dir, sep='\t')

