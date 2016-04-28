__author__ = 'dp'
import pandas
import re
import operator
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle
import matplotlib
import copy
import HTSeq

norm_left_len = 200.
norm_right_len = 200.
norm_cds_len = 1000.
utr_len = 800
verbose = False


class flocs:
    """Feature locations for a transcript.
    """
    def __init__(self, transcript_id, chrom, txpt_left, txpt_right,
                 cds_left, cds_right, strand, seq, exons_dict,
                 gene_name, exon_borders_in_seq, no_left_utr, no_right_utr):
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.txpt_left = txpt_left
        self.txpt_right = txpt_right
        self.cds_left = cds_left
        steps_backward = 0
        init_cds_left = self.cds_left
        self.start_pos_in_seq = self.cds_left - self.txpt_left + 1
        while (seq[self.start_pos_in_seq:self.start_pos_in_seq+3] != 'ATG'):
            if cds_left < 1: break
            self.cds_left -= 1
            steps_backward += 1
            if steps_backward > 4: break
            self.start_pos_in_seq = self.cds_left - self.txpt_left + 1
        if (seq[self.start_pos_in_seq:self.start_pos_in_seq+3] != 'ATG'):
            self.cds_left = init_cds_left
        self.cds_right = cds_right
        self.strand = strand
        #if self.strand == '+':
        #    self.cds_left -= 1
        self.seq = seq  # Upper case.
        self.exons_dict = exons_dict
        self.gene_name = gene_name
        self.exon_borders_in_seq = exon_borders_in_seq
        # self.fbe_locs = []
        self.motif_locs = []
        self.polyA_locs = []
        self.fivep_seq = ''
        self.threep_seq = ''
        self.peak_locs = []
        self.peak_locs_genomic_coord = []
        if verbose: print "type cds_left %s txpt_left %s" % (type(self.cds_left), type(self.txpt_left))
        self.left_len = float(self.cds_left - self.txpt_left)
        self.right_len = float(self.txpt_right - self.cds_right)
        self.cds_len = len(self.seq) - self.left_len - self.right_len
        self.stop_pos_in_seq = len(seq) - (self.txpt_right - self.cds_right) # Right border of trinucleotide seq.
#s        if self.gene_name == 'mrpl-15':
        if verbose: print "flocs.__init__(): %s, len(seq) %i, txpt_right %f, cds_right %f, stop_pos_in_seq %f" % (
            self.gene_name, len(seq), self.txpt_right, self.cds_right, self.stop_pos_in_seq)
        self.start_pos_in_seq = self.cds_left - self.txpt_left + 1 # Left border of trinucleotide seq.
        self.norm_neg1 = []
        self.norm_neg2 = []
        self.norm_neg3 = []
        self.neg1_locs = []
        self.neg2_locs = []
        self.neg3_locs = []
        self.neg1_vs_polyA = []
        self.neg2_vs_polyA = []
        self.neg3_vs_polyA = []

    def check_coherency(self):
        li = """{name}: start[{a0}-{a1}]: {start} stop[{b0}-{b1}]: {stop} \
given left utr: {given_l_utr} given right utr: {given_r_utr} \
seq len: {seqlen}""".format(
            name=self.gene_name,
            start=self.seq[self.start_pos_in_seq:self.start_pos_in_seq+3],
            a0=self.start_pos_in_seq, a1=self.start_pos_in_seq+3,
            stop=self.seq[self.stop_pos_in_seq:self.stop_pos_in_seq+3],
            b0=self.stop_pos_in_seq, b1=self.stop_pos_in_seq+3,
            given_l_utr=self.given_left_utr, given_r_utr=self.given_right_utr,
            seqlen=len(self.seq))
        if self.seq[self.start_pos_in_seq:self.start_pos_in_seq+3] != 'ATG':
            print "Wrong start:\n" + li
            print "flocs.__init__(): %s, len(seq) %i, txpt_right %f, cds_right %f, stop_pos_in_seq %f" % (
                self.gene_name, len(self.seq), self.txpt_right, self.cds_right, self.stop_pos_in_seq)
            print str(self)
            print "Left UTR:\n" + self.seq[0:self.start_pos_in_seq]
            print "CDS:\n" + self.seq[self.start_pos_in_seq:self.stop_pos_in_seq]
            print self.seq
            print '***'
        if self.seq[self.stop_pos_in_seq:self.stop_pos_in_seq+3] not in [
            'TAA', 'TAG', 'TGA']:
            print "Wrong stop:\n" + li
            print str(self)
            print "Left UTR:\n" + self.seq[0:self.start_pos_in_seq]
            print "CDS:\n" + self.seq[self.start_pos_in_seq:self.stop_pos_in_seq]
            print "Right UTR:\n" + self.seq[self.stop_pos_in_seq:]
            print self.seq
            print '***'

    def __str__(self):
        li = """
name: {name}
\tcds_start: {s} cds_end: {e} txpt_start: {ts} txpt_end: {te}
\tstrand: {strand}
\tmotif_locs: {fbes}
\tpolyA_locs: {polya}
\tpeak_locs: {pks}
\tseq_len: {seqlen}
\tend of seq: {endseq}
""".format(
    name=self.gene_name,
    s=self.cds_left, e=self.cds_right,
    ts=self.txpt_left,
    te=self.txpt_right,
    strand=self.strand,
    fbes=str(self.motif_locs), polya=str(self.polyA_locs),
    pks=str(self.peak_locs),
    seqlen=len(self.seq),
    endseq=self.seq[-10:])
        if hasattr(self, 'norm_motif'):
            li += "\tnorm fbes: %s\n" % str(self.norm_motif)
        if hasattr(self, 'norm_peaks'):
            li += "\tnorm peaks: %s\n" % str(self.norm_peaks)
        if hasattr(self, 'norm_polyA'):
            li += "\tnorm polyA: %s\n" % str(self.norm_polyA)
        return li

    def find_features(self, peaks=False):
        print 'finding features'
        self.find_motif()
        print self.motif_no_star
        print self.motif_locs

        self.find_controls()
        self.find_polyA()
        if peaks:
            self.find_peaks(peaks)

    def find_motif(self):
        #self.seq = self.seq.lower()
        # p = re.compile('TGT[^C][^G][^G]AT')
        self.motif_no_star = self.motif
        if '\*'  in self.motif:
            self.motif_no_star = re.sub('\*', '', self.motif_no_star)
        if '.' in self.motif_no_star:
            self.motif_no_star = re.sub('\.', '', self.motif_no_star)
        self.motif_no_star = '(' + self.motif_no_star + ')'
        p = re.compile(self.motif, re.IGNORECASE)
        print self.seq

        for m in p.finditer(self.seq):
            self.motif_locs.append(
                (m.start(), m.end())
                )

    def find_controls(self):
        #self.seq = self.seq.lower()
        p1 = re.compile('ATGATG')  # Use A\wGTT\w\wT
        p2 = re.compile('A\wGTT\w\wT')  #AGA\w\wT\wA seems to be 3'UTR enriched like the FBE. Unclear why.
        p3 = re.compile('G\wAT\w\wTTA')
        for m in p1.finditer(self.seq):
            self.neg1_locs.append(
                (m.start(), m.end())
                )
        for m in p2.finditer(self.seq):
            self.neg2_locs.append(
                (m.start(), m.end())
                )
        for m in p3.finditer(self.seq):
            self.neg3_locs.append(
                (m.start(), m.end())
                )
        if verbose: print "find_controls(): %s %s %s" % (str(self.neg1_locs), str(self.neg2_locs), str(self.neg3_locs))

    def find_polyA(self):
        p = re.compile('AATAAA')
        for m in p.finditer(self.seq):
            if m.start() > self.stop_pos_in_seq:
                self.polyA_locs.append(
                    (m.start(), m.end())
                    )
        if len(self.polyA_locs) > 0:
            self.last_polyA = sorted(self.polyA_locs, key=lambda x: x[0])[-1]
        else: self.last_polyA = False

    def set_primary_and_secondary(self):
        if len(self.peak_locs) > 1:
            highest_peak = sorted(self.peak_locs, key=lambda x: float(x[2]))[-1]

    def get_highest_peak(self):
        return sorted(self.peak_locs_genomic_coord, key=lambda x: float(x[2]))[-1]

    def get_distances_between_peaks(self):
        self.intrapeak_distances = []
        if len(self.peak_locs) < 2:
            self.intrapeak_distances = [] #np.array([])
        else:
            for index_a, peak_a in enumerate(self.peak_locs):
                for index_b, peak_b in enumerate(self.peak_locs):
                    if index_a >= index_b: continue
                    self.intrapeak_distances.append(
                        self.dist_between_two_peaks(peak_a, peak_b)
                    )
            #self.intrapeak_distances = np.array(self.intrapeak_distances)

    def dist_between_two_peaks(self, a, b):
        if a[0] <= b[0] <= a[1]:
            return 0
        if a[0] <= b[1] <= a[1]:
            return 0
        if b[0] >= a[1]:
            return b[0] - a[1]
        elif b[1] <= a[0]:
            return a[0] - b[1]

    def add_raw_reads_to_utr(self, ga, chr_len):
        #utr_left = self.cds_right
        #utr_right = self.txpt_right
        # Need chrom information.
        self.utr_arr = []
        if self.txpt_right - self.cds_right < 2:
            for pos in range(0, self.txpt_right - self.cds_right + 1, 1):
                self.utr_arr = [0]
            return
        if self.strand == '-':
            txpt_left = chr_len[self.chrom] - self.txpt_right + 1
            txpt_right = chr_len[self.chrom] - self.txpt_left
            cds_left = chr_len[self.chrom] - self.cds_right + 1
            cds_right = chr_len[self.chrom] - self.cds_left
        else:
            txpt_left = self.txpt_left
            txpt_right = self.txpt_right
            cds_left = self.cds_left
            cds_right = self.cds_right
        iv = HTSeq.GenomicInterval(
            self.chrom, cds_right,
            txpt_right, self.strand)
        self.utr_ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
        # if len(ga[iv].steps()) == 0:
        #    for pos in range(0,self.txpt_right - self.cds_right + 1,1):
        #        self.utr_arr.append(0)
        #    return
        if txpt_right - cds_right < 2:
            self.utr_arr = [0]
            return
        for _iv, score in ga[iv].steps():
            self.utr_ga[_iv] = score
            left_in_utr = _iv.start - cds_right
            right_in_utr = _iv.end - cds_right
            for pos in range(left_in_utr, right_in_utr, 1):
                self.utr_arr.append(score)

    def add_raw_fbes_to_all_peak_regions(self, seqs):
        self.all_peaks_fbes_in_region = []
        for peak in self.peak_locs_genomic_coord:  # Left, right, height.
            center_of_peak = int(float(peak[0] + peak[1])/2.)
            left = center_of_peak - 1000
            right = center_of_peak + 1000
            gen_seq_in_region = seqs[self.chrom][left:right]
            if self.strand == '-':
                gen_seq_in_region = rc(gen_seq_in_region)

    def add_raw_reads_to_all_peak_regions(self, ga):
        self.highest_peak_raw = False
        self.highest_peak_highest_marked = []
        self.secondary_peaks_raw = []
        self.secondary_peaks_highest_marked = []
        self.region_around_highest_raw = False
        self.all_peaks_raw = []
        self.all_peaks_highest_marked = []
        self.region_around_all_peaks_raw = []
        if len(self.peak_locs) == 0: return
        if len(self.peak_locs) == 1:
            highest_peak_height = self.peak_locs[0][2]
        if len(self.peak_locs) > 1:
            highest_peak_height = sorted(self.peak_locs, key=lambda x: float(x[2]))[-1][2]
        for peak in self.peak_locs_genomic_coord:  # left, right, height
            peak_raw = self.add_raw_reads_to_a_peak_region(peak, ga)
            highest_range = self.get_highest_peak() # left, right, height
            marks = self.add_raw_reads_to_a_peak_region(peak, ga, range_to_mark=highest_range)
            if peak_raw == [0]: continue
            (left_raw, right_raw) = self.add_raw_reads_outside_peak_region(peak, ga)
            if left_raw == [0]: continue
            self.all_peaks_raw.append(peak_raw)
            self.all_peaks_highest_marked.append(marks)
            self.region_around_all_peaks_raw.append([left_raw, right_raw])
            if peak[2] == highest_peak_height:
                self.highest_peak_raw = peak_raw
                self.highest_peak_highest_marked.append(marks)
                self.region_around_highest_raw = [left_raw, right_raw]
            if peak[2] < highest_peak_height:
                self.secondary_peaks_raw.append(peak_raw)
                self.secondary_peaks_highest_marked.append(marks)

    def add_raw_reads_to_a_peak_region(self, peak, ga, range_to_mark=None):
        peak_center = int(float(peak[1] + peak[0])/2.)
        #print "Adding raw reads to a peak in %s" % self.gene_name
        left_raw = []
        right_raw = []
        if peak_center <= 1000:
            left_border = 1
            to_pad = peak_center - 1000 - 1
            left_raw = [0] * to_pad
        else:
            left_border = peak_center - 1000
        right_border = peak_center + 1000
        if peak_center - left_border < 2: return [0]
        if range_to_mark is not None:
            marks = []
            left_marks_border = max(range_to_mark[0], left_border)
            right_marks_border = min(range_to_mark[1], right_border)
            # Case 1: the range to mark does not overlap the peak range.
            if not (
                        (left_border <= range_to_mark[0] <= right_border) or (
                        left_border <= range_to_mark[1] <= right_border)
                    ):
                marks += [0] * (int(right_border) - int(left_border))
            else:

                marks += [0] * (left_marks_border - left_border)
                marks += [1] * (right_marks_border - left_marks_border)
                marks += [0] * (right_border - right_marks_border)
            return marks
            nope = '''
            # Case 2: the range to mark is within the peak range.
            elif (
                        (left_border <= range_to_mark[0] <= right_border) and (
                        left_border <= range_to_mark[1] <= right_border)
                    ):
                marks += [0] * (left_marks_border - left_border)
                marks += [1] * (right_marks_border - left_marks_border)
                marks += [0] * (right_border - right_marks_border)
            # Case 3: only the right range overlaps.
            elif (
                    (left_border > range_to_mark[0]) and (
                        left_border <= range_to_mark[1] <= right_border)
                    ):

            )
            if int(left_border) < range_to_mark[0] < int(right_border):
                marks += [0] * (range_to_mark[0] - int(left_border))
                marks += [1] * (min(range_to_mark[1], int(right_border)) - range_to_mark[0])
            if range_to_mark[1] > range_to_mark[0]:
                marks += [1] * (range_to_mark[1] - range_to_mark[0])
            if peak_center + 1000 > range_to_mark[1]:
                marks += [0] * (peak_center + 1000 - range_to_mark[1])
            return marks'''
        left_iv = HTSeq.GenomicInterval(
            self.chrom, int(left_border), peak_center, self.strand)
        right_iv = HTSeq.GenomicInterval(
            self.chrom, peak_center, peak_center + 1000, self.strand)
        for iv, score in ga[left_iv].steps():
            left_raw += [score] * (iv.end - iv.start)
        for iv, score in ga[right_iv].steps():
            right_raw += [score] * (iv.end - iv.start)
        peak_raw = left_raw + right_raw
        return peak_raw

    def add_raw_reads_outside_peak_region(self, peak, ga):
        left_raw = []
        right_raw = []
        if peak[0] <= 1000:
            left_border = 1
            to_pad = peak[0] - 1000 - 1
            left_raw = [0] * to_pad
        else:
            left_border = peak[0] - 1000
        if peak[0] - left_border < 2: return [0], [0]
        left_iv = HTSeq.GenomicInterval(
            self.chrom, left_border, peak[0], self.strand)
        right_iv = HTSeq.GenomicInterval(
            self.chrom, peak[1], peak[1] + 1000, self.strand)
        for iv, score in ga[left_iv].steps():
            left_raw += [score] * (iv.end - iv.start)
        for iv, score in ga[right_iv].steps():
            right_raw += [score] * (iv.end - iv.start)
        return [left_raw, right_raw]

    def find_peaks(self, peaks_l, chr_lens):
        'peaks : list of dict objects.'
        if verbose: print "flocs.find_peaks(): Looking for %s." % self.gene_name
        if self.gene_name not in peaks_l:
            return False
        # Peaks on - strand have coordinates flipped, as just like the floc object.
        for row in peaks_l[self.gene_name]:
            if (self.txpt_left <= row['left'] <= self.txpt_right) or (
                    self.txpt_left <= row['right'] <= self.txpt_right):
                if row['right'] > self.txpt_right:
                    peak_right = self.txpt_right
                else: peak_right = row['right']
                if row['left'] < self.txpt_left:
                    peak_left = self.txpt_left
                else: peak_left =row['left']
                # Find which exon the left border is in.
                self.peak_exon = 0
                left_found = 0
                right_found = 0
                for exon_num in self.exons_dict:
                    try:
                        exon = self.exons_dict[exon_num]
                    except:
                        if verbose: print 'flocs.find_peaks(): Did not find %s in %s. Exons_dict = %s' % (
                            str(exon_num), str(self.exons_dict), str(self.exons_dict))
                    if (exon[1] <= peak_left <= exon[2]):
                        # Where in the mature sequence is this exon?
                        (exon_start_pos_in_seq, end) = self.exon_borders_in_seq[exon_num]
                        pos_in_seq_left = exon_start_pos_in_seq + (peak_left - exon[1])
                        left_found = 1
                        if verbose: print "flocs.find_peaks(): Found peak left border %i in exon %i (%s)." % (peak_left, exon_num, str(exon))
                        if verbose: print "exon_start_pos_in_seq, exon[1] = %i, %i" % (exon_start_pos_in_seq, exon[1])
                    if (exon[1] <= peak_right <= exon[2]):
                        (exon_start_pos_in_seq, end) = self.exon_borders_in_seq[exon_num]
                        pos_in_seq_right = exon_start_pos_in_seq + (peak_right - exon[1])
                        right_found = 1
                        if verbose: print "flocs.find_peaks(): Found peak right border %i in exon %i (%s)." % (peak_right, exon_num, str(exon))
                try:
                    if left_found and right_found:
                        self.peak_locs.append(
                            (pos_in_seq_left, pos_in_seq_right, float(row['height'])))
                        # - strand peaks need their coordinates flipped to compare
                        # with the flocs object. However, afterwards, it is more
                        # useful to use them with true coordinates (for the sake of
                        # obtaining raw signal around the peak).
                        if self.strand == '-':
                            new_peak_left = chr_lens[self.chrom] - peak_right + 1
                            peak_right = chr_lens[self.chrom] - peak_left
                            peak_left = new_peak_left
                        self.peak_locs_genomic_coord.append(
                            (peak_left, peak_right, float(row['height']))
                        )
                        if verbose:
                            print "flocs.find_peaks(): Appended %s to self.peak_locs." % str((
                                pos_in_seq_left, pos_in_seq_right, row['height']))
                            print "flocs.find_peaks(): chrom %s, added %s to self.peak_locs_genomic_coord" % (
                                self.chrom, str((peak_left, peak_right, float(row['height'])))
                            )
                except:
                    #if verbose:
                    print "Peak not found in exon: %i. txpt: %i %i" % (
                        peak_left, int(self.txpt_left), int(self.txpt_right))
                    #if verbose:
                    print "Exons: %s" % str(self.exons_dict)
                    #if verbose:
                    print "Exon borders in seq dict: %s" % str(self.exon_borders_in_seq)
            else:  # Peak not within transcript.
                pass

    def normalize_features(self):
        if verbose: print "flocs.normalize_features() called for %s" % self.transcript_id
        self.norm_motif = []
        self.norm_polyA = []
        self.norm_peaks = []
        if verbose: print "All data in this object: %s" % str(self.__dict__)
        for fbe in self.motif_locs:
            self.norm_motif.append(
                (self.get_norm_pos(fbe[0],test='fbe'), self.get_norm_pos(fbe[1],test='fbe'))
                )
        for (norm_set, motif_set) in [(self.norm_neg1, self.neg1_locs),
                                      (self.norm_neg2, self.neg2_locs), (self.norm_neg3, self.neg3_locs)]:
            for motif in motif_set:
                norm_set.append(
                    (self.get_norm_pos(motif[0],test='motif'), self.get_norm_pos(motif[1],test='motif'))
                    )
        for peak in self.peak_locs:
            if verbose: print "flocs.normalize_features(): peak %s" % str(peak)
            self.norm_peaks.append(
                (self.get_norm_pos(peak[0], test='peak'), self.get_norm_pos(peak[1],test='peak'), peak[2])  # Third element is height.
                )
            if verbose: print "flocs.normalize_features(): Normalized peak %s to %s. CDS borders %i, %i." % (
                str(peak), str(self.norm_peaks[-1]), self.start_pos_in_seq, self.stop_pos_in_seq)
        for polyA in self.polyA_locs:
            self.norm_polyA.append(
                (self.get_norm_pos(polyA[0],test='polyA'), self.get_norm_pos(polyA[1],test='polyA'))
                )
        if verbose: print "normalize_features(): %s %s %s" % (str(self.norm_neg1), str(self.norm_neg2), str(self.norm_neg3))

    def get_norm_pos(self, pos, test=''):
        if pos < self.start_pos_in_seq:
            norml = int(float(pos) * norm_left_len/float(self.start_pos_in_seq))
        if float(self.start_pos_in_seq) <= pos <= float(self.stop_pos_in_seq):
            pos_in_cds = pos - self.start_pos_in_seq
            norml = int(
                norm_left_len + int(float(pos_in_cds) * norm_cds_len/float(self.stop_pos_in_seq - self.start_pos_in_seq))
                )
        if pos > float(self.stop_pos_in_seq):
            pos_in_right = pos - self.stop_pos_in_seq
            norml = int(
                norm_left_len + norm_cds_len + int(float(pos_in_right) * norm_right_len/float(len(self.seq) - self.stop_pos_in_seq))
                )
        if not hasattr(self, 'motif'):
            print 'Motif not set. Set floc.motif before calling get_norm_pos.'
        # if motif is None:
        #     motif = 'TGT\w\w\wAT'
        try:
            if norml > 1400:
                if verbose: print "\ntype %s" % test
                if verbose: print "\tpos %i left_len %f cds_len %f right_len %f sum of lens %f seq len %i" % (
                    pos, self.left_len, self.cds_len, self.right_len,
                    self.left_len + self.cds_len + self.right_len,
                    len(self.seq))
                if verbose: print "\tnorml = %i" % norml
                p = re.compile(motif)
                for m in p.finditer(self.seq):
                    if verbose: print m.span()
                    #print self.get_norm_pos(m.span()[0])
            return norml
        except:
            if verbose: print "Error in flocs.get_norm_pos(): pos=%s" % str(pos)
            return 0

    def get_peak_pos_relative_to_polyA(self):
        if len(self.polyA_locs) == 0: return
        #self.norm_peak_v_polyA = [0.,0.]
        if len(self.peak_locs) > 0:
            #print "Found a peak."
            highest_peak = sorted(self.peak_locs, key=lambda x: float(x[2]))[-1]
            self.norm_peak_v_polyA = self.get_norm_pos_of_feature(highest_peak)
            if self.norm_peak_v_polyA:
                self.norm_peak_v_polyA.append(highest_peak[2])  # Add height.

    def get_motif_pos_relative_to_polyA(self):
        if len(self.motif_ocs) == 0: return
        self.norm_motif_v_polyA = []
        if len(self.motif_locs) > 0:
            for fbe in self.motif_locs:
                if type(fbe) != tuple:
                    print "fbe not tup: %s" % str(self.motif_locs)
                    continue
                _norm_p = self.get_norm_pos_of_feature(fbe)
                if _norm_p:
                    self.norm_motif_v_polyA.append(_norm_p)

    def get_motif_pos_relative_to_polyA(self):
        for (vs_polyA, norm_set) in [
            (self.neg1_vs_polyA, self.norm_neg1), (self.neg2_vs_polyA, self.norm_neg2),
            (self.neg3_vs_polyA, self.norm_neg3)]:
            if len(norm_set) == 0: return
            if len(norm_set) > 0:
                self.vs_polyA = []
                for motif in norm_set:
                    if type(motif) != tuple:
                        print "motif not tup: %s" % str(norm_set)
                        continue
                    _norm_p = self.get_norm_pos_of_feature(motif)
                    if _norm_p:
                        vs_polyA.append(_norm_p)

    def get_norm_pos_of_feature(self, feat):
        norm_CDS_to_polyA = 200.
        norm_polyA_to_end = 100.
        norm_feat = [0., 0.]
        if len(self.polyA_locs) < 1: return
        last_polyA = self.last_polyA
        if feat[0] < self.stop_pos_in_seq:
            return
        dist_peak_to_polyA = last_polyA[0] - feat[0]
        # Normalize to UTR length.
        utr_len = float(len(self.seq) - self.stop_pos_in_seq)
        dist_CDS_to_polyA = last_polyA[0] - self.stop_pos_in_seq
        dist_polyA_to_end = len(self.seq) - last_polyA[0]
        norm_CDS_to_polyA_factor = norm_CDS_to_polyA/float(dist_CDS_to_polyA)
        norm_polyA_to_end_factor = norm_polyA_to_end/float(dist_polyA_to_end)
        if feat[0] < last_polyA[0]:
            norm_feat[0] = (feat[0] - self.stop_pos_in_seq) * norm_CDS_to_polyA_factor
        else:
            norm_feat[0] = norm_CDS_to_polyA + (feat[0] - last_polyA[0]) * norm_polyA_to_end_factor
        if feat[1] < last_polyA[0]:
            norm_feat[1] = (feat[1] - self.stop_pos_in_seq) * norm_CDS_to_polyA_factor
        else:
            norm_feat[1] = norm_CDS_to_polyA + (feat[1] - last_polyA[0]) * norm_polyA_to_end_factor
        return norm_feat
