import HTSeq
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import pickle
import pandas
import re
import os
#import pysam
import scipy as sp
#import statsmodels.api as sm
#from rpy2.robjects.packages import importr
#import rpy2
#MASS = importr('MASS')
#rstats = importr('stats')


def check_overlap(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda x: x.left)
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher.left <= lower.right:
                print "Overlapping peak ranges..\n%s\n%s\n*" % (str(lower), str(higher))
        

def get_distance_pandas(_gene, apeak):
    """Used to assign peaks to a gene.
    """
    if _gene['4'] < apeak.left:
        return apeak.left - _gene['4']
    if apeak.right < _gene['3']:
        return apeak.right - _gene['3']
    return 0


def find_nearby_genes(sub, apeak):
    """Assigns a gene for a given peak.
    """
    cutoff = 1e3
    genes_for_peak = []
    for index, _gene in sub.iterrows():
        dist = get_distance_pandas(_gene, apeak)
        if abs(dist) < cutoff:
            genes_for_peak.append(_gene)
    return genes_for_peak


def resolve_tie(ties):
    biotypes = []
    for tup in ties:
        _gene = tup[1]
        m = re.search('gene_biotype "([^"]+)"', _gene['8'])
        if m is not None:
            biotypes.append((tup[0], tup[1], m.group(1)))
        else:
            biotypes.append((tup[0], tup[1], 'Unknown'))
    # If there is a non-coding RNA, assign to that.
    non_coding = []
    for tup in biotypes:
        if tup[2] != 'protein_coding':
            non_coding.append(tup)
    if len(list(non_coding)) == 1:
        return (non_coding[0][0], non_coding[0][1])
    if len(list(non_coding)) > 1:
        # Multiple non-coding RNAs. Pick randomly.
        return (non_coding[0][0], non_coding[0][1])
    if len(list(non_coding)) == 0:
        # No non-coding RNA. Pick randomly.
        return (ties[0][0], ties[0][1])


# Find the closet gene for each peak, and resolve ties.   
def assign_to_gene(merged_peaks, chrm, strand, gtf_df, verbose=False):
    if chrm not in merged_peaks:
        return False
    if strand not in merged_peaks[chrm]:
        return False
    sub = gtf_df[(gtf_df['0']==chrm) & (gtf_df['6']==strand)]
    sub = sub[sub['2']=='transcript']
    #sub_list_left = list(sub['3'])  # List of left positions.
    #sub_list_right = list(sub['4'])  # List of left positions.
    for index, apeak in enumerate(merged_peaks[chrm][strand]):
        if not index % 100:
            print "Assigning gene %i/%i for chrm/strand %s/%s" % (
                index, len(merged_peaks[chrm][strand]), chrm, strand)
        #index_left = bisect.bisect_left(sub_list_left, apeak.left - 1e3)
        #index_right = bisect.bisect_left(sub_list_right, apeak.left + 1e3)
        asub = sub[(abs(sub['3'] - apeak.left) < 1e5
                    ) | (abs(sub['4'] - apeak.left) < 1e5)]
        apeak.genes_for_peak = find_nearby_genes(asub, apeak)
        if verbose:
            print "assign_to_gene(): .genes_for_peak="
            for a_gene in apeak.genes_for_peak:
                print "\n%s\n" % str(a_gene)
        closest = (1e4, None)
        ties = []
        for _gene in apeak.genes_for_peak:
            dist = get_distance_pandas(_gene, apeak)
            if verbose:
                print "\n***_gene %s: dist: %i***\n" % (_gene.gene_name, dist) 
            if abs(dist) < abs(closest[0]):
                ties = []
                closest = (dist, _gene)
                if verbose:
                    print "Set closest to %s, dist %i" % (_gene.gene_name, dist)
            elif abs(dist) == abs(closest[0]):
                ties.append((dist, dict(_gene)))
                ties.append((dist, dict(closest[1])))
        if len(ties) > 0:
            for tup in ties:
                dist = tup[0]
                _gene = tup[1]
                #dist = get_distance_pandas(_gene, apeak)
                if dist == closest[0]:
                    closest = resolve_tie(ties)
                    break
        apeak.gene_dist = closest[0]
        apeak.gene = closest[1]
        del apeak.genes_for_peak


def create_gtf_with_names_file(gtf_noheader_filename):
    outf = open('gtf_with_names_column.txt', 'w')
    with open(gtf_noheader_filename, 'r') as f:
        outf.write('\t'.join([str(x) for x in range(0,9)]) + '\tgene_name\ttranscript_id\ttranscript_name\texon_number\n')
        for li in f:
            s = li.rstrip('\n').split('\t')
            s = s + [get_name(li), get_txpt_name(li), get_txpt_id(li), get_exon_num(li)]
            outf.write("\t".join(s) + '\n')
    outf.close()
    gtf_df = pandas.read_csv('gtf_with_names_column.txt', sep='\t', header=0)
    return gtf_df


def get_name(line):
    #line = str(line)
    m = re.search('gene_name "([^"]+)"', line)
    if m is not None:
        return m.group(1)
    else:
        return 'NA'


def get_txpt_id(line):
    #line = str(line)
    m = re.search('transcript_id "([^"]+)"', line)
    if m is not None:
        return m.group(1)
    else:
        return 'NA'


def get_txpt_name(line):
    #line = str(line)
    m = re.search('transcript_name "([^"]+)"', line)
    if m is not None:
        return m.group(1)
    else:
        return 'NA'


def get_exon_num(line):
    #line = str(line)
    m = re.search('exon_number "([^"]+)"', line)
    if m is not None:
        return m.group(1)
    else:
        return 'NA'


def complement(s):
    #basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s


# For each exon, get the 5' read ends.
def read_region(bam_filename, iv):
    bamfile = pysam.AlignmentFile(bam_filename, "rb")
    s = bamfile.fetch(iv[0], max(0, int(iv[1])), iv[2])
    reads = list()
    ga = HTSeq.GenomicArray([iv[0]], stranded=True)
    for r in s:
        if(iv[3]=="+" and not r.is_reverse):
            r_pos = HTSeq.GenomicPosition(
                 iv[0], r.reference_start, iv[3])
            ga[r_pos] += 1
        if (iv[3]=="-" and r.is_reverse):
           r_pos = HTSeq.GenomicPosition(
                 iv[0], r.reference_end-1, iv[3])
           ga[r_pos] += 1
    bamfile.close()
    return ga

def total_reads_in_bin(ga, iv):
    total_reads = 0
    for iv, value in list(ga[iv].steps()):
        total_reads += value
    return total_reads

def split_into_bins(ga, region):
    bin_size = 50
    bins_list = []
    for bin_left in range(region[1], region[2], bin_size):
        bin_right = bin_left + bin_size
        iv = HTSeq.GenomicInterval(region[0], bin_left, bin_right, region[3])
        bins_list.append(total_reads_in_bin(ga, iv))
    return bins_list

def add_local_signal(peak_objs, bamfiles):
    print "Adding local signal..."
    genes = {}
    width = 500
    for peak in peak_objs:
        for bamfile in bamfiles:  # Values are filenames.
            # Determine local signal.
            region = [peak.chrm,
                      max(0, peak.left-width),
                      peak.right+width,
                      peak.strand] # Chr, start, end, strand.
            reads = read_region(str(bamfiles[bamfile]), region)
            bins_list = split_into_bins(reads, region)
            peak.local[bamfile] = bins_list


def add_gene_signal(peak_objs, gtf_df, bamfiles):
    print "Adding gene signal..."
    genes = {}
    for peak in peak_objs:
        if peak.gene is None or 'transcript_id' not in peak.gene:
            continue
        txpt_id = peak.gene['transcript_id']
        exons_df = gtf_df[(gtf_df['transcript_id']==txpt_id) & (gtf_df['2']=='exon')]
        exons = zip(exons_df['0'], exons_df['3'], exons_df['4'], exons_df['6'])
        for bamfile in bamfiles:  # Values are filenames.
            if txpt_id in genes and bamfile in genes[txpt_id]:
                peak.exons[bamfile] = genes[txpt_id][bamfile]
            else:
                # Determine gene signal.
                bins_list = []
                if txpt_id not in genes:
                    genes[txpt_id] = {}
                peak.exons[bamfile] = []
                genes[txpt_id][bamfile] = []
                for exon in exons:  # Each element is an interval tuple.
                    reads = read_region(bamfiles[bamfile], exon)
                    bins_list = split_into_bins(reads, exon)
                    genes[txpt_id][bamfile].extend(bins_list)
                    peak.exons[bamfile].extend(bins_list)


def do_statistics(peak_objs, bamfiles):
    bamfile_sizes = {}
    norm_factors = {}
    for bamfile in bamfiles:
        # One-liner from https://www.biostars.org/p/1890/:
        bamfile_sizes[bamfile] = reduce(
            lambda x, y: x + y,
            [eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamfiles[bamfile])])
    norm_factors['clip'] = 1.0
    norm_factors['neg_ip'] = float(bamfile_sizes['clip'])/float(
        bamfile_sizes['neg_ip'])
    norm_factors['rna_seq'] = float(bamfile_sizes['clip'])/float(
        bamfile_sizes['rna_seq'])
    for peak in peak_objs:
        peak.find_max_bin()
        peak.local_poisson(bamfiles, norm_factors)
        peak.local_norm(bamfiles, norm_factors)
        peak.gene_poisson(bamfiles, norm_factors)
        peak.gene_norm(bamfiles, norm_factors)
