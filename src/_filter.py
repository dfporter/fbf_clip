"""

A second version of the config file settings would be the config.ini format:

[library]

top: ./
lib: /groups/Kimble/Common/fog_iCLIP/calls/lib/
gtf_raw: %(lib)s/Caenorhabditis_elegans.WBcel235.78.noheader.gtf
fasta: %(lib)s/c_elegans.WS235.genomic.fa
fai: %(lib)s/c_elegans.WS235.genomic.fa.fai
chr_sizes: %(lib)s/chr_sizes.txt
gtf_one_txpt_per_gene: %(lib)s/gtf_one_txpt_per_gene.txt
gtf: %(lib)s/gtf_with_names_column.txt

bedgraphs_folder: %(top)s/../bedgraphs_unnorm/combined_unnorm/

# Used by subpeaks.py
bedgraph_exp_plus:  %(bedgraphs_folder)s/exp_+.wig
bedgraph_exp_minus: %(bedgraphs_folder)s/exp_-.wig
read_beds:  %(top)s/../bed_collapsed/combined/
#
figs: %(top)s/figs/
#
control_bed1: control_n2.bed
#control_bed2: n2_oo_lane1_rt16.bed
#control_bed3: n2_oo_lane1_rt3.bed
exp_bed1: fbf_rep_1.bed
exp_bed2: fbf1_oo_lane2_rt6.bed
exp_bed3: fbf1_oo_lane2_rt9.bed
#
clusters_dir: %(top)s/clusters/
permutation_peaks_dir: %(top)s/permutation_peaks/
"""
import re
import os
import sys
import pandas
import argparse
import HTSeq
import collections
import numpy as np
import config
import os


def get_val(ga, iv):
    return np.max(np.fromiter(ga[iv], dtype=np.float))


def add_heights_to_peak_file(peak_fname, bed_ga_dict):
    if type(peak_fname) == type(''):
        try:
            peaks = pandas.read_csv(peak_fname, sep='\t')
        except:
            print "Could not open {a}. Empty file?".format(a=peak_fname)
            return pandas.DataFrame()
    else:
        peaks = peak_fname
        pass # Assume dataframe. Lazy.
    print "Height columns will appear as %s \n=\n%s.\n " % (
        bed_ga_dict.keys(), [x + '_reads' for x in bed_ga_dict])
    ivs = zip(peaks['chrm'].tolist(), peaks['left'].tolist(),
        peaks['right'].tolist(), peaks['strand'].tolist())
    for bedgraph_name in bed_ga_dict:
        peaks[bedgraph_name] =  [
            get_val(bed_ga_dict[bedgraph_name], HTSeq.GenomicInterval(*iv)
                    ) for iv in ivs]
        peaks[bedgraph_name + '_reads'] = peaks[bedgraph_name]
    return peaks


def load_bedgraph_file(fname, add_strand=True):
    fname = fname.rstrip('.bed')
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    if add_strand:
        plus_file = fname.partition('.wig')[0] + '_+.wig'
        add_strand_to_ga_from_bedgraph_file(plus_file, ga, '+')
        minus_file = fname.partition('.wig')[0] + '_-.wig'
        add_strand_to_ga_from_bedgraph_file(minus_file, ga, '-')
    else:
        if re.match('.*\+.*', os.path.basename(fname)) is not None:
            add_strand_to_ga_from_bedgraph_file(fname, ga, '+')
        elif re.match('.*-.*',  os.path.basename(fname)) is not None:
            add_strand_to_ga_from_bedgraph_file(fname, ga, '-')
    return ga


def add_strand_to_ga_from_bedgraph_file(fname, ga, strand):
    with open(fname, 'r') as f:
        next(f)
        for line in f:
            #if line[0] != 'M': continue
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), strand)] = float(s[3])
    return ga


def get_bed_size(fname):
    return float(len(open(fname).readlines()))


def add_read_columns(args, config, bedfiles=None):
    non_control = {}
    if bedfiles is None:
        bedfiles = {'control': "{a}/{b}.wig".format(
            a=os.path.dirname(config['clip_replicate'][0]),
            b=os.path.basename(config['neg_ip_filename']).rstrip('.bed'))}
        for x in config['clip_replicate']:
            name = os.path.basename(x).rstrip('.bed').rstrip('.wig')
            bedfiles[name] = x
            non_control[name] = x
    ga_d = {}
    print "Loading bedgraphs from:"
    for k in bedfiles: print "{a}: {b}".format(a=k, b=bedfiles[k])
    for name in bedfiles:
        ga_d[name] = load_bedgraph_file(bedfiles[name])
    peaks = add_heights_to_peak_file(args.peaks_fname, ga_d)
    peaks = normalize_height_columns(bedfiles, config, peaks)
    return peaks


def normalize_height_columns(bedfiles, config, peaks):
    size = {}
    for bedgraph in bedfiles:
        if re.search('control', bedgraph):
            size[bedgraph] = get_bed_size(config['neg_ip_filename'])
            col_name = 'depth_control_' + bedgraph
            peaks[col_name] = [
                1e6 * x/size[bedgraph] for x in peaks[bedgraph].tolist()]
            continue
        print config['bed_dir']
        print bedgraph
        print os.path.basename(bedgraph).rstrip('.wig')
        bed_file = config['bed_dir'] +'/' + \
                   os.path.basename(bedgraph).rstrip('.wig') + '.bed'
        size[bedgraph] = get_bed_size(bed_file)
        col_name = 'depth_exp_' + bedgraph
        peaks[col_name] = [
            1e6 * x/size[bedgraph] for x in peaks[bedgraph].tolist()]
    return peaks


def add_sum_and_ratio_columns(peaks):
    exp_cols = [col for col in peaks.columns if (
        re.search('\Aexp_', col) and re.search('_reads$', col))]
    control_cols = [col for col in peaks.columns if (
        re.search('\Acontrol_', col) and re.search('_reads$', col))]
    print "Using columns %s and %s for ratio calculation." % (
        exp_cols, control_cols)
    peaks['ratio'] = -1
    peaks['exp'] = 0
    peaks['control'] = 0
    for i, row in peaks.iterrows():
        exp = sum([
            peaks.loc[i, col] for col in exp_cols])/float(len(exp_cols))
        control = sum([
            peaks.loc[i, col] for col in control_cols])/float(len(control_cols))
        ratio = exp/max([1., float(control)])
        peaks.loc[i, 'exp'] = exp
        peaks.loc[i, 'control'] = control
        peaks.loc[i, 'ratio'] = ratio
        #print '{a}/{b} = {c}'.format(
        #    a=exp, b=control, c=ratio)
        #peaks.loc[i, 'ratio'] = -100
    return peaks


def get_bedgraph_to_bed_size(lib):
    bedgraph_basenames = dict([(os.path.basename(x).partition('.wig')[0], x) \
        for x in lib['clip_replicate']])
    bed_basenames = dict([(os.path.basename(x).partition('.bed')[0], x) \
        for x in lib['clip_replicate_bed']])
    print bedgraph_basenames
    print bed_basenames
    bedgraph_to_bed = dict([(bedgraph_basenames[x], bed_basenames[x]) for \
        x in bedgraph_basenames])
    bedgraph_negative = "{a}/{b}.wig".format(a=lib['bedgraphs_folder'], 
        b=os.path.basename(lib['neg_ip_filename'].partition('.bed')[0]))
    print bedgraph_to_bed
    bed_size = dict([(x, sum([1 for x in open(bedgraph_to_bed[x]).readlines()])) \
        for x in bedgraph_to_bed])
    bed_size[bedgraph_negative] = sum([1 for x in \
        open(lib['neg_ip_filename']).readlines()])
    print bed_size
    return bed_size

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
Loads an unnormalized bedgraph file and determines the maximum coverage\
 in each peak region. Bed files are loaded and used to normalize the \
read depths to per billion. Ratio cutoffs are applied to this ratio.''')
    parser.add_argument('-p', '--peaks_fname',
                        help='Filename of peaks file.')
    parser.add_argument('-c', '--config')
    parser.add_argument('-r', '--ratio_cutoff',
                        help='Enrichment ratio cutoff.')
    parser.add_argument('-o', '--output', default='output.txt')
    args = parser.parse_args()
    args.ratio_cutoff = float(args.ratio_cutoff)
    lib = config.config(filepath=args.config)
    header = open(args.peaks_fname).readline()
    if (re.search('depth_exp', header) is not None) and (
        re.search('depth_control', header) is not None):
        print "Peaks file appears to already have read count columns.\
 Overwrite them?"
        answer = raw_input('(Y/N?) >')
        print answer
        if answer[0].upper() == 'Y':
            peaks = add_read_columns(args, lib)
            peaks = add_sum_and_ratio_columns(peaks)
            peaks.to_csv(args.peaks_fname, sep='\t', index=False)
        else:
            print "Using the existing columns then..."
    else:
        peaks = add_read_columns(args, lib)
        peaks = add_sum_and_ratio_columns(peaks)
        peaks.to_csv(args.peaks_fname, sep='\t', index=False)
    peaks = pandas.read_csv(args.peaks_fname, sep='\t')
    peaks = peaks[peaks['ratio']>=args.ratio_cutoff]
    if not os.path.exists(os.path.dirname(args.output)):
        if len(os.path.dirname(args.output)):
            os.system('mkdir ' + os.path.dirname(args.output))
    peaks.to_csv(args.output, sep='\t', index=False)

