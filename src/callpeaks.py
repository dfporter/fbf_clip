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
import glob
#import pysam
import scipy as sp
#import statsmodels.api as sm
from peakslib import *
#from statsmodels.sandbox.stats.multicomp import multipletests


def find_peaks(coverage, chrm='I', strand='+'):
    #if chrm != 'I':
    #    return []
    for iv in coverage[chrm][strand].steps():
        pass
    last_value = iv[0].start
    min_read_cutoff = 10
    peak_pos = []
    width = 1e3
    for start in range(0, last_value, int(width)):
    #for start in range(15050000, 15060000, int(width)):
        end = start + width
        window = HTSeq.GenomicInterval(chrm, start, end, strand)
        coverage_array = np.fromiter(coverage[window], dtype='i')
        if max(coverage_array) < min_read_cutoff:
            continue
        xs = np.arange(100, 200, 100)
        peakind = signal.find_peaks_cwt(coverage_array, xs, noise_perc=20)
        for a_pos in peakind:
            peak_pos.append(start + a_pos)
        #if start > 1e5:
        #    break
    return peak_pos

# From: http://stackoverflow.com/questions/17118350/how-to-find-nearest-value-that-is-greater-in-numpy-array
def argfind(array, predicate):
    for i in xrange(array.shape[0]):
        if predicate(array[i]):
            return i
    return False


def find_nearest_below(array, value):
    return argfind(array, lambda x: x <= value)


def find_borders(peak_pos, coverage, chrm, strand):
    # What are the borders of each peak?
    peak_objs = []
    for a_peak_pos in peak_pos:
        width = 1e3
        left = max(0, a_peak_pos - width)
        right = a_peak_pos + width
        window = HTSeq.GenomicInterval(chrm, left, right, strand)
        peak_center = HTSeq.GenomicPosition(chrm, a_peak_pos, strand)
        if coverage[peak_center] < 5:
            continue
        wincvg = np.fromiter(coverage[window], dtype='i')
        cutoff = max(wincvg)*.1
        rv = wincvg[::-1]
        left_border = find_nearest_below(rv[1000:], cutoff)
        right_border = find_nearest_below(wincvg[1000:], cutoff)
        left_border = a_peak_pos - left_border
        right_border = a_peak_pos + right_border
        #print "find_borders(): peak_pos %s:%i" % (chrm, peak)
        #print "peak value %s" % str(coverage[peak_center])
        #print "%s:%i-%i" % (chrm,left_border,right_border)
        
        if left_border == right_border:
            left_border -= 10
            right_border += 10
            #print "Error, left=right border: %s" % str(peak_obj.__dict__)
            #window = HTSeq.GenomicInterval(chrm, left_border-10, right_border+10, strand)
            #wincvg = np.fromiter(coverage[window], dtype='i')
            #peak_obj.height = int(max(wincvg))
        peak_obj = peak(chrm, left_border, right_border, strand)
        window = HTSeq.GenomicInterval(chrm, left_border, right_border, strand)
        wincvg = np.fromiter(coverage[window], dtype='i')
        peak_obj.height = int(max(wincvg))
        if peak_obj.height == 'na':
            print "Set a peak height to 'na': %s..." % str(peak_obj.__dict__)
        peak_objs.append(peak_obj)
    return peak_objs

        
# Merge overlapping peaks.
# From: http://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
def merge_overlapping_on_chrm_and_strand(intervals, coverage):
    """Merge in a given chrom and strand.
    """
    sorted_by_lower_bound = sorted(intervals, key=lambda x: x.left)
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher.left <= lower.right:
                upper_bound = int(max(lower.right, higher.right))
                new_peak = peak(lower.chrm, lower.left, upper_bound, lower.strand)
                window = HTSeq.GenomicInterval(lower.chrm, lower.left, upper_bound, lower.strand)
                wincvg = np.fromiter(coverage[window], dtype='i')
                new_peak.height = int(max(wincvg))
                merged[-1] = new_peak  # replace by merged interval
            else:
                merged.append(higher)
    return merged

def merge_overlapping(intervals):
    by_chrm = {}
    merged = {}
    for apeak in intervals:
        by_chrm.setdefault(apeak.chrm, {})
        merged.setdefault(apeak.chrm, {})
        by_chrm[apeak.chrm].setdefault(apeak.strand, [])
        merged[apeak.chrm].setdefault(apeak.strand, [])
        by_chrm[apeak.chrm][apeak.strand].append(apeak)
    for chrm in by_chrm:
        for strand in by_chrm[chrm]:
            merged[chrm][strand] = merge_overlapping_on_chrm_and_strand(by_chrm[chrm][strand])
            # Check.
            check_overlap(merged[chrm][strand])
    return merged


def any_have_na(peak_objs):
    has_na = False
    for apeak in peak_objs:
        if apeak.height == 'na':
            print "Has na height: %s" % str(apeak.__dict__)
    if not has_na:
        print "No peaks have na height."


class peak:
    def __init__(self, chrm, left, right, strand):
        self.chrm = str(chrm)
        self.left = int(left)
        self.right = int(right)
        self.strand = str(strand)
        self.iv = HTSeq.GenomicInterval(chrm, left, right, strand)
        self.local = {}
        self.exons = {}
        self.pvalues = {}
        
    def __repr__(self):
        #return "__dict__:%s" % str(self.__dict__)
        if hasattr(self, 'height'):
            return "%s:%i-%i (height=%i)" % (self.chrm, self.left, self.right, self.height)
        return "%s:%i-%i" % (self.chrm, self.left, self.right)

    def find_max_bin(self):
        """The peak bin is the maximum bin in the center 40% of bins.
        40% * 1kbp = 400 bp = 8 (50 bp) bins.
        """
        x = self.local['clip']
        #print "x: %s\n length: %i" % (str(x), len(x))
        #print "left boundary index: %i" % int(len(x)*0.2)
        #print "right bountdary index: %i" % int(len(x)*.2 + int(len(x)*.5))
        midrange = x[int(len(x)*0.2):int(len(x)*.2 + int(len(x)*.5))]
        #print "Searching midrange %s for max_bin." % str(midrange)
        self.max_bin = max(midrange)
        #print "Set max bin as %s." % str(self.max_bin)

    def local_poisson(self, bamfiles, norm_factors):
        mean_local = {}
        for bamfile in bamfiles:
            mean_local[bamfile] = float(sum(self.local[bamfile]))/float(
                len(self.local[bamfile]))
            mean_local[bamfile] *= norm_factors[bamfile]
            self.pvalues['%s_local_poisson' % bamfile] = sp.stats.poisson.sf(
                self.max_bin, mean_local[bamfile])
            if np.isnan(self.pvalues['%s_local_poisson' % bamfile]):
                if mean_local[bamfile] == 0.0 and self.max_bin > 0:
                    self.pvalues['%s_local_poisson' % bamfile] = 0
                if mean_local[bamfile] == 0.0 and self.max_bin == 0:
                    self.pvalues['%s_local_poisson' % bamfile] = 1.0
            # Bonferroni to make p value the significance for local region.
            self.pvalues['%s_local_poisson' % bamfile] *= len(self.local[bamfile])
            self.pvalues['%s_local_poisson' % bamfile] = min(1., self.pvalues['%s_local_poisson' % bamfile])
            if np.isnan(self.pvalues['%s_local_poisson' % bamfile]):
                print "Error in local_poisson(): %s" % str(self.__dict__)

    def local_norm(self, bamfiles, norm_factors):
        mean_local = {}
        for bamfile in bamfiles:
            norm_list = [norm_factors[bamfile] * x for x in self.local[bamfile]]
            self.local_norm_mu, self.local_norm_std = sp.stats.norm.fit(norm_list)
            self.pvalues['%s_local_norm' % bamfile] = 1 - sp.stats.norm.cdf(
                self.max_bin, self.local_norm_mu, self.local_norm_std)
            if np.isnan(self.pvalues['%s_local_norm' % bamfile]):
                if self.local_norm_mu == 0 and self.max_bin > 0:
                    self.pvalues['%s_local_norm' % bamfile] = 0
                elif self.local_norm_mu == 0 and self.max_bin == 0:
                    self.pvalues['%s_local_norm' % bamfile] = 1.0
            # Bonferroni to make p value the significance for the gene.
            self.pvalues['%s_local_norm' % bamfile] *= len(self.local[bamfile])
            self.pvalues['%s_local_norm' % bamfile] = min(1., self.pvalues['%s_local_norm' % bamfile])
                
    def gene_poisson(self, bamfiles, norm_factors):
        mean_gene = {}
        for bamfile in bamfiles:
            if not hasattr(self, 'exons') or bamfile not in self.exons:
                if not hasattr(self, 'exons'):
                    print "No exons for %s." % str(self.__dict__)
                elif bamfile not in self.exons:
                    print "No bamfile %s in self.exons %s" % (bamfile, str(self.exons))
                self.pvalues['%s_gene_poisson' % bamfile] = 'na'
                continue
            mean_gene[bamfile] = float(sum(self.exons[bamfile]))/float(
                len(self.exons[bamfile]))
            mean_gene[bamfile] *= norm_factors[bamfile]
            self.pvalues['%s_gene_poisson' % bamfile] = sp.stats.poisson.sf(
                self.max_bin, mean_gene[bamfile])
            if np.isnan(self.pvalues['%s_gene_poisson' % bamfile]):
                if mean_gene[bamfile] == 0.0 and self.max_bin > 0:
                    self.pvalues['%s_gene_poisson' % bamfile] = 0
                elif mean_gene[bamfile] == 0.0 and self.max_bin == 0:
                    self.pvalues['%s_gene_poisson' % bamfile] = 1.0
            # Bonferroni to make p value the significance for the gene.
            self.pvalues['%s_gene_poisson' % bamfile] *= len(self.exons[bamfile])
            self.pvalues['%s_gene_poisson' % bamfile] = min(1., self.pvalues['%s_gene_poisson' % bamfile])

    def gene_norm(self, bamfiles, norm_factors):
        for bamfile in bamfiles:
            if not hasattr(self, 'exons') or bamfile not in self.exons:
                self.gene_norm_mu, self.gene_norm_std = ('na', 'na')
                self.pvalues['%s_gene_norm' % bamfile] = 'na'
                continue
            norm_list = [norm_factors[bamfile] * x for x in self.exons[bamfile]]
            self.gene_norm_mu, self.gene_norm_std = sp.stats.norm.fit(norm_list)
            self.pvalues['%s_gene_norm' % bamfile] = 1 - sp.stats.norm.cdf(
                self.max_bin, self.gene_norm_mu, self.gene_norm_std)
            if np.isnan(self.pvalues['%s_gene_norm' % bamfile]):
                if self.gene_norm_mu == 0 and self.max_bin > 0:
                    self.pvalues['%s_gene_norm' % bamfile] = 0
                elif self.gene_norm_mu == 0 and self.max_bin == 0:
                    self.pvalues['%s_gene_norm' % bamfile] = 1.0
            # Bonferroni to make p value the significance for the gene.
            self.pvalues['%s_gene_norm' % bamfile] *= len(self.exons[bamfile])
            self.pvalues['%s_gene_norm' % bamfile] = min(1., self.pvalues['%s_gene_norm' % bamfile])

def check_bins(peak_objs, bamfiles):
    print "Checking %i peaks..." % len(peak_objs)
    for apeak in peak_objs:
        if not hasattr(apeak, 'local'):
            print "\nNo .local attribute for %s." % str(apeak.__dict__)
            continue
        for bamfile in bamfiles:
            if bamfile not in apeak.local:
                print "\nNo bamfile %s in local for %s." % (bamfile, str(apeak.__dict__))
        if not hasattr(apeak, 'exons'):
            print "\nNo .exons attribute for %s." % str(apeak.__dict__)
            continue
        for bamfile in bamfiles:
            if bamfile not in apeak.exons:
                print "\nNo bamfile %s in exons for %s." % (bamfile, str(apeak.__dict__))

def check_for_nan(peak_objs, bamfiles):
    skip = 1
    indx = 0
    for p in peak_objs:
        row = [p.chrm, p.left, p.right, p.strand, p.height, p.max_bin]
        if not hasattr(p, 'gene') or p.gene is None:
            continue
        if len(p.pvalues) != 9:
            print "Less than 9 values for peak %s" % str(p.__dict__)
            if indx > skip:
                return p
            else:
                indx += 1
            continue
        row += [str(p.gene[key]) for key in dict(p.gene).keys()]
        pvalues = [str(p.pvalues[key]) for key in dict(p.pvalues).keys()]
        for x in p.pvalues.values():
            if np.isnan(x):
                print "Found nan value in %s" % str(p.__dict__)
                mean_gene = {}
                for bamfile in bamfiles:
                    print "bamfile: %s" % bamfile
                    mean_gene[bamfile] = float(sum(p.exons[bamfile]))/float(len(p.exons[bamfile]))
                    if np.isnan(p.pvalues['%s_local_poisson' % bamfile]):
                        if mean_local[bamfile] == 0.0 and p.max_bin > 0:
                            p.pvalues['%s_local_poisson' % bamfile] = 0
                            print "mean local == 0 and max bin > 0"
                        elif mean_local[bamfile] == 0.0 and p.max_bin == 0:
                            p.pvalues['%s_local_poisson' % bamfile] = 1.0
                            print "mean local == 0 and max bin == 0"
                        else:
                            print "mean local != 0 or max bin < 0"
    print "Checked for nan in %i peaks." % len(peak_objs)
#p = check_for_nan(peak_objs_by_chrm['I']['+'], bamfiles)
#for chrm in peak_objs_by_chrm:
#    for strand in peak_objs_by_chrm[chrm]:
        #do_statistics(peak_objs_by_chrm[chrm][strand], bamfiles)
#        check_for_nan(peak_objs_by_chrm[chrm][strand], bamfiles)

def call_peaks_on_replicate(clip_bam_filename, config, load_data=False):
    loaded_data = False
    if load_data:
        try:
            with open('peak_objs_by_chrm_%s.p' % os.path.basename(clip_bam_filename), 'rb') as f:
                peak_objs_by_chrm = pickle.load(f)
            peak_table = convert_peak_objs_to_table(peak_objs_by_chrm)
            loaded_data = True
            print "Loaded peaks called previously and saved as file %s." % str(
                'peak_objs_by_chrm_%s.p' % os.path.basename(clip_bam_filename))
        except:
            print "Failed to load data for %s. Regenerating..." % clip_bam_filename
            loaded_data = False
    if not loaded_data:
        peak_objs_by_chrm = call_peaks_from_bam(clip_bam_filename, config)
        peak_table = convert_peak_objs_to_table(peak_objs_by_chrm)
        with open('peak_objs_by_chrm_%s.p' % os.path.basename(clip_bam_filename), 'wb') as f:
            pickle.dump(peak_objs_by_chrm, f)
    bamfiles = {'clip': clip_bam_filename,
                'rna_seq': config['rna_seq_filename'],  # .bam filename.
                'neg_ip': config['neg_ip_filename']}  # .bam filename.
    fdr_correction(bamfiles, peak_table)
    # Write peaks and evaluate for each of the seven different null hypothesis.
    evaluate_hypothesis(peak_table, clip_bam_filename)
    return peak_table

def call_peaks_from_bam(clip_bam_filename, config):
    gtf_filename = config['gtf_filename'] #"/home/dp/Desktop/celegans_genome/wormbase_ws235/Caenorhabditis_elegans.WBcel235.78.gtf"
    gtf_noheader_filename = config['gtf_filename_noheader'] #"/home/dp/Desktop/celegans_genome/wormbase_ws235/Caenorhabditis_elegans.WBcel235.78.noheader.gtf" 
    gtffile = HTSeq.GFF_Reader(gtf_filename)
    #clip_bam_filename = "/home/dp/Desktop/bams/celegans/run813_fbf_aacc_20mapq.bam"
    bamfiles = {'clip': clip_bam_filename,
                'rna_seq': config['rna_seq_filename'],
                'neg_ip': config['neg_ip_filename']}
    clip_bamfile = HTSeq.BAM_Reader(clip_bam_filename)
    coverage = HTSeq.GenomicArray("auto", stranded=True, typecode='i')
    #gtf_df = pandas.read_csv(gtf_noheader_filename, sep='\t', header=None)
    print "Reading alignments from bamfile..."
    for aln in clip_bamfile:  # Very slow.
        if aln.aligned:
            coverage[aln.iv] += 1
    print "Creating gtf file and dataframe..."
    gtf_df = create_gtf_with_names_file(gtf_noheader_filename)
    print "Calling peaks..."
    peaks_by_chrm = {}
    peak_objs_by_chrm = {}
    for chrm in dict(gtf_df['0'].value_counts()).keys():
        peaks_by_chrm[chrm] = {}
        peak_objs_by_chrm[chrm] = {}
        for strand in ['+', '-']:
            peaks_by_chrm[chrm][strand] = find_peaks(coverage, chrm=chrm, strand=strand)
            peak_objs_by_chrm[chrm][strand] = find_borders(
                peaks_by_chrm[chrm][strand], coverage, chrm, strand)
            peak_objs_by_chrm[chrm][strand] = merge_overlapping_on_chrm_and_strand(
                peak_objs_by_chrm[chrm][strand], coverage)
            assign_to_gene(peak_objs_by_chrm, chrm, strand, gtf_df)
            add_local_signal(peak_objs_by_chrm[chrm][strand], bamfiles)
            add_gene_signal(peak_objs_by_chrm[chrm][strand], gtf_df, bamfiles)
            do_statistics(peak_objs_by_chrm[chrm][strand], bamfiles)
            any_have_na(peak_objs_by_chrm[chrm][strand])
    return peak_objs_by_chrm

def convert_peak_objs_to_table(peak_objs_by_chrm):
    # Convert to a table and write.
    peak_list = []
    output_cols = ['chrm', 'left', 'right', 'strand',
                   'height', 'max_bin']  # Simple columns.
    output_cols += dict(peak_objs_by_chrm['I']['+'][0].gene).keys()
    output_cols += dict(peak_objs_by_chrm['I']['+'][0].pvalues).keys()
    for chrm in peak_objs_by_chrm:
        for strand in peak_objs_by_chrm[chrm]:
            for p in peak_objs_by_chrm[chrm][strand]:
                row = [p.chrm, p.left, p.right, p.strand, p.height, p.max_bin]
                if not hasattr(p, 'gene') or p.gene is None:
                    continue
                row += [str(p.gene[key]) for key in dict(p.gene).keys()]
                row += [str(p.pvalues[key]) for key in dict(p.pvalues).keys()]
                #else:
                #    row += ['na' for x in range(1,14)]
                #    row += ['na' for x in range(1,10)]
                peak_list.append(row)
    peak_table = pandas.DataFrame(peak_list, columns=output_cols)
    #no_gene = peak_table[peak_table['transcript_id']=='na']
    peak_table = peak_table[peak_table['transcript_id']!='na']
    peak_table.to_csv('peak_table.txt', sep='\t')
    return peak_table


def fdr_correction(bamfiles, peak_table):
    # Need FDRs.
    fdrs = {}
    for bamfile in bamfiles:
        # Local Poisson (only use CLIP).
        fdrs["%s_local_poisson" % bamfile] = multipletests(peak_table["%s_local_poisson" % bamfile].astype(float),
                                      alpha=0.05, method='fdr_bh')
        peak_table["%s_local_poisson_rej" % bamfile] = pandas.Series(
            fdrs["%s_local_poisson" % bamfile][0], index=peak_table.index)
        peak_table["%s_local_poisson_cor" % bamfile] = pandas.Series(
            fdrs["%s_local_poisson" % bamfile][1], index=peak_table.index)
        # Gene Poisson (only use CLIP).
        fdrs["%s_gene_poisson" % bamfile] = multipletests(peak_table["%s_gene_poisson" % bamfile].astype(float),
                                      alpha=0.05, method='fdr_bh')
        peak_table["%s_gene_norm_rej" % bamfile] = pandas.Series(
            fdrs["%s_gene_poisson" % bamfile][0], index=peak_table.index)
        peak_table["%s_gene_poisson_cor" % bamfile] = pandas.Series(
            fdrs["%s_gene_poisson" % bamfile][1], index=peak_table.index)
        # Local normal.
        fdrs["%s_local_norm" % bamfile] = multipletests(peak_table["%s_local_norm" % bamfile].astype(float),
                                      alpha=0.05, method='fdr_bh')
        peak_table["%s_local_norm_rej" % bamfile] = pandas.Series(
            fdrs["%s_local_norm" % bamfile][0], index=peak_table.index)
        peak_table["%s_local_norm_cor" % bamfile] = pandas.Series(
            fdrs["%s_local_norm" % bamfile][1], index=peak_table.index)
        # Gene normal.
        fdrs["%s_gene_norm" % bamfile] = multipletests(peak_table["%s_gene_norm" % bamfile].astype(float),
                                      alpha=0.05, method='fdr_bh')
        peak_table["%s_gene_norm_rej" % bamfile] = pandas.Series(
            fdrs["%s_gene_norm" % bamfile][0], index=peak_table.index)
        peak_table["%s_gene_norm_cor" % bamfile] = pandas.Series(
            fdrs["%s_gene_norm" % bamfile][1], index=peak_table.index)

def evaluate_hypothesis(peak_table, clip_bam_filename):
    replicate_dir = os.path.basename(clip_bam_filename).rstrip('.bam')
    os.system('mkdir %s' % replicate_dir)
    # Null hypothesis 2: Signal is not a local peak relative to itself.
    sub = peak_table[peak_table['clip_local_poisson_cor']<0.05]
    sub.to_csv('%s/null_hyp_2.txt' % replicate_dir, sep='\t')
    # Null hypothesis 3: Signal is not a peak relative to CLIP signal in the gene.
    sub = peak_table[peak_table['clip_gene_poisson_cor']<0.05]
    sub.to_csv('%s/null_hyp_3.txt' % replicate_dir, sep='\t')
    # Null hypothesis 4: Signal is not enriched relative to the negative in the local region.
    sub = peak_table[peak_table['neg_ip_local_norm_cor']<0.05]
    sub.to_csv('%s/null_hyp_4.txt' % replicate_dir, sep='\t')
    # Null hypothesis 5: Signal is not enriched relative to the negative for the gene.
    sub = peak_table[peak_table['neg_ip_gene_norm_cor']<0.05]
    sub.to_csv('%s/null_hyp_5.txt' % replicate_dir, sep='\t')
    # Null hypothesis 6: Signal is not enriched relative to RNA-seq locally.
    sub = peak_table[peak_table['rna_seq_local_norm_cor']<0.05]
    sub.to_csv('%s/null_hyp_6.txt' % replicate_dir, sep='\t')
    # Null hypothesis 7: Signal is not enriched relative to RNA-seq for the gene.
    sub = peak_table[peak_table['rna_seq_gene_norm_cor']<0.05]
    sub.to_csv('%s/null_hyp_7.txt' % replicate_dir, sep='\t')

def read_config(filename):
    """Expect:
experiment_name\tname  # Optional.
clip_replicate\tbamfilename1
clip_replicate\tbamfilename2
clip_replicate\tbamfilename3
gtf_filename\tgtf_filename
gtf_filename_noheader\tfilename
rna_seq_filename\tfilename
neg_up_filename\tfilename
    """
    config = {}
    with open(filename, 'r') as f:
        for li in f:
            li = li.partition('#')[0]  # Skip comments.
            s = li.rstrip('\n').split('\t')
            config.setdefault(s[0], [])
            try: config[s[0]].append(s[1])
            except: print "Error parsing line %s in config file. Wrong length?" % li
    for key in config:
        if len(config[key]) == 1:
            config[key] = config[key][0]
    if 'experiment_name' not in config:
        config['experiment_name'] = os.path.basename(filename)
    return config

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

def score_metrics(dir_name, config):
    li = ""
    for filename in glob.glob(dir_name + '/*'):
        print filename
        li += score_metric(filename)
    print li
    with open('score_metrics_%s.txt' % config['experiment_name'], 'w') as f:
        f.write(li)
    return li

def score_metric(filename, label="", given_peaks=False, peaks=False):
    if not label:
        label = os.path.basename(filename)
    if not given_peaks:
        peaks = pandas.DataFrame.from_csv(filename, sep='\t')
    get_sequences(peaks)
    score_binding_site(peaks)
    run_dreme(peaks, label)
    positives = score_positives(peaks)
    return write_metrics(peaks, positives, label)

def run_dreme(combined, label):
    if not os.path.exists('combined_fasta'):
        os.system('mkdir combined_fasta')
    fasta_filename = 'combined_fasta/%s.fa' % label
    as_fasta = ""
    for index, peak_row in combined.iterrows():
        start = combined.loc[index, 'left']
        end = combined.loc[index, 'right']
        chrm = combined.loc[index, 'chrm']
        name = combined.loc[index, 'gene_name']
        name += ":%s:%i:%i" % (chrm, start, end)
        as_fasta += ">%s\n%s\n" % (name, combined.loc[index, 'seq'])
    with open(fasta_filename, 'w') as f:
        f.write(as_fasta)
    out_dir = 'dreme_results/%s' % label
    os.system('/home/dp/meme/bin/dreme -e 0.000000001 -p %s -norc -maxk 10 -oc %s' % (fasta_filename, out_dir))

def get_sequences(combined):
    fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for index, peak_row in combined.iterrows():
        start = combined.loc[index, 'left']
        end = combined.loc[index, 'right']
        chrm = combined.loc[index, 'chrm']
        seq = sequences[chrm][start:end]
        if combined.loc[index, 'strand'] == '-':
            seq = rc(seq)
        combined.loc[index, 'seq'] = seq

def write_metrics(peaks, positives, label):
    li = """
Dataset: {label}
Number of peaks: {df_size}
Number of genes: {n_genes}
With FBE: {with_fbe}, {fbe_perc}%
Without FBE: {without_fbe}
Positive controls: {observed}/{expected}
Missing positives: {missing}
""".format(label=label,
    df_size=len(peaks), n_genes=len(list(set(peaks['gene_name']))),
           with_fbe=len(peaks[peaks['has_fbe']==1]),
           fbe_perc= float(100 * len(peaks[peaks['has_fbe']==1])/len(peaks)),
           without_fbe=len(peaks[peaks['has_fbe']==0]),
           observed=positives['number of observed positives'],
           expected=positives['expected'],
           missing=positives['missing positives'])
    return li

def score_binding_site(peaks):
    for index, peak_row in peaks.iterrows():
        #print 'searching %s' % peaks.loc[index, 'seq']
        if re.search('tgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'has_fbe'] = 1
        else:
            peaks.loc[index, 'has_fbe'] = 0

def score_positives(peaks):
    known_pos = set(['gld-1', 'htp-1', 'htp-2', 'mpk-1', 'him-3',
                         'fbf-1', 'lip-1', 'syp-2', 'fbf-2', 'fog-1',
                         'fem-3', 'syp-3', 'gld-3', 'fog-3', 
                         'egl-4'])
    obs_genes = set(peaks['gene_name'])
    obs_pos = known_pos & obs_genes
    missing_pos = known_pos - obs_genes
    obs_pos_n = len(list(obs_pos))
    missing_pos_n = len(list(missing_pos))
    return {'observed positives': obs_pos, 'number of observed positives': obs_pos_n,
            'missing positives': missing_pos, 'number of missing positives': missing_pos_n,
            'expected': len(list(known_pos))}

def combine_replicates(config):
    combined = {}
    hyp = {}
    rep_dir_names = [
        "./%s/" % os.path.splitext(os.path.basename(rep_bam_path))[0] for rep_bam_path in config['clip_replicate']]
    for hyp_num in range(2, 8):
        hyp[hyp_num] = {}
        for rep_name in rep_dir_names:
            filename = '%s/null_hyp_%i.txt' % (rep_name, hyp_num)
            hyp[hyp_num][rep_name] = get_peaks(filename)
        combined[hyp_num] = combine_peaks(hyp[hyp_num])
        write_combined(combined[hyp_num], str(hyp_num), config)
    return combined

def get_peaks(filename):
    peaks = pandas.DataFrame.from_csv(filename, sep='\t')
    return peaks

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
        if len(overlapping_peaks) > 1:
            consensus = consensus_peak(overlapping_peaks)
            tup = (consensus['chrm'], consensus['left'], consensus['strand'])
            if tup not in combined:
                combined[tup] = consensus
        #if len(combined) > 10:
        #    return combined
    return combined

def consensus_peak(peaks):
    rows = peaks[peaks['height']==max(peaks['height'])]
    return rows.iloc[0].to_dict()


#def overlapping(iv, comb_df):
#    sub = comb_df[(comb_df['chrm']==iv[0]) & (comb_df['strand']==iv[3])]
#    overlap = sub[(sub['left']<=iv[2]) & (iv[1]<=sub['right'])]
#    return overlap


def overlapping(iv, comb_df):
    overlap = comb_df[(comb_df['chrm']==iv[0]) & (comb_df['strand']==iv[3]) & (comb_df['left']<=iv[2]) & (iv[1]<=comb_df['right'])]
    return overlap


def write_combined(combined, label, config):
    combined = pandas.DataFrame(combined.values())
    if not os.path.exists('combined'):
        os.system('mkdir combined_%s/' % config['experiment_name'])
    with open('combined_%s/%s' % (config['experiment_name'], label), 'w') as f:
        combined.to_csv(f, sep='\t')


if __name__ == '__main__':
    config = read_config('config_fbf1')
    for clip_replicate in config['clip_replicate']:
        call_peaks_on_replicate(clip_replicate, config, load_data=True)
    combine_replicates(config)
    score_metrics('combined_%s/' % config['experiment_name'], config)

