"""
     txpt objects, bedgraph coverage    peaks file, sequence data
                \                         \
subpeaks_ui.py ---- subpeaks.run()-------------process_file()

process_file()
|
|_________analyze_by_applying_cwt_to_the_peak_regions()
|
|_________analyze_by_applying_cwt_to_all_reads_in_gene()

analyze_by_applying_cwt_to_the_peak_regions()
|
|_for each peak, define a superpeak, split by CWT
   \_peak.subpeak_stats()
      \_Create info dict and equivalent _stats object
      |  \_pass all peaks to stats_on_all_superpeaks
      |     \_stats_on_all_superpeaks() combines _stats objects and writes result (averages...)
      |         \_pass lists of superpeaks to image creation functions: pulsar, heatmaps
      \___Write fastas.

analyze_by_applying_cwt_to_all_reads_in_gene()
|
|_for each gene target, create sbGene object of the gene
   \_call CWT and pass save rows of unnormalized signal around subpeaks as array
      \_create subpeak_rows object for combined info on highest peaks and secondaries.
         |_pass normalized rows from the subpeak_rows to the pulsar graphing function
         |_call the graphing functions of the subpeak_rows objects for average peak graphs

"""

import os
import sys

import HTSeq
import pandas
import build_image_objects_for_heatmaps
import output_heatmap_figures
import logging
import datetime
import time
import numpy as np
import matplotlib
from scipy import signal
import pulsar
import re

from scipy.signal import argrelextrema
logger = logging.getLogger()
import matplotlib.pyplot as plt
from find_local_extrema import find_local_extrema

import edge_of_txpt
class subpeak:

    def __init__(self, iv):
        self.iv = iv


class superpeak:
    """
    find_subpeaks(coverage)-->subpeak_arr2d [[position, height]]
    Outside superpeak, convert subpeak_arr2d objects to genomic coordinates,
    combine all subpeak_arr2d objects (for all superpeak objects) for a gene, and sort.

    """
    def __init__(self, iv):
        self.iv = iv

    def find_subpeaks(self, coverage):
        width = 0
        expanded_iv = HTSeq.GenomicInterval(
            self.iv.chrom,
            max([self.iv.start - width, 1]),
            self.iv.end + width, self.iv.strand,
        )
        self.peak_arr = np.fromiter(coverage[expanded_iv], dtype='i')
        self.sub_peaks = self.do_cwt(self.peak_arr)
        self.subpeak_ivs_relative = [HTSeq.GenomicPosition(
            self.iv.chrom, pos, self.iv.strand,
            ) for pos in self.sub_peaks]
        self.subpeak_arr2d = [[pos, self.peak_arr[pos]] for pos in self.sub_peaks]
        self.subpeak_arr2d = sorted(self.subpeak_arr2d, key=lambda x: x[1])
        if len(self.subpeak_arr2d) > 0:
            self.highest_subpeak = self.subpeak_arr2d[-1]
            self.subpeak_stats()
        else:
            self.highest_subpeak = None

    def give_rows_of_subpeaks_for_image(self, coverage, chr_lens):
        width = 1e3
        self.subpeak_rows = []
        if not hasattr(self, 'unnormalized_subpeak_rows'):
            self.unnormalized_subpeak_rows = self.give_unnormalized_rows_of_subpeaks_for_image(
                coverage, chr_lens)
        for unnormalized_subpeak_row in self.unnormalized_subpeak_rows:
            a_subpeak_row = normalize_row_to_highest(unnormalized_subpeak_row)
            a_subpeak_row = np.clip(a_subpeak_row, 1e-9, 1)
            self.subpeak_rows.append(a_subpeak_row)
        return self.subpeak_rows

    def get_seqs_from_ivs(self, sequences):
        if sequences is None: return []
        self.seqs = []
        if len(self.subpeak_arr2d) < 1: return []
        for _iv in self.get_subpeak_ivs_given_width(width=35):
            if _iv.end > len(sequences[_iv.chrom]):
                end = len(sequences[_iv.chrom])
            else: end = _iv.end
            start = max([_iv.start, 0])
            seq = sequences[_iv.chrom][start:end]
            if _iv.strand == '-':
                seq = rc(seq)
            self.seqs.append(seq)  # Highest peak last.

    def get_subpeak_ivs_given_width(self, width=50):
        subpeak_region_ivs = []
        self.subpeak_arr2d = sorted(self.subpeak_arr2d, key=lambda x: x[1])  # Should be already sorted.
        for pos in self.subpeak_arr2d:
            center_genomic = self.iv.start + pos[0]
            subpeak_region_ivs.append(HTSeq.GenomicInterval(
                self.iv.chrom, max([1, center_genomic - width]),
                center_genomic + width, self.iv.strand))
        return subpeak_region_ivs

    def get_subpeak_region_ivs(self):
        width = 1e3
        self.subpeak_ivs = []
        if len(self.subpeak_arr2d) < 1: return []
        self.subpeak_arr2d = sorted(self.subpeak_arr2d, key=lambda x: x[1])  # Should be already sorted.
        for pos in self.subpeak_arr2d:
            center_genomic = self.iv.start + pos[0]
            self.subpeak_ivs.append(HTSeq.GenomicInterval(
                self.iv.chrom, max([1, center_genomic - width]),
                center_genomic + width, self.iv.strand))
        return self.subpeak_ivs

    def give_unnormalized_rows_of_subpeaks_for_image(self, coverage, chr_lens):
        self.unnormalized_subpeak_rows = []
        for subpeak_iv in self.get_subpeak_region_ivs():
            unnormalized_subpeak_row = np.nan_to_num(np.fromiter(coverage[subpeak_iv], dtype='i'))
            unnormalized_subpeak_row = np.clip(unnormalized_subpeak_row, 1e-9, 1e9)
            self.unnormalized_subpeak_rows.append(unnormalized_subpeak_row)
            self.unnormalized_subpeak_rows = sorted(
                self.unnormalized_subpeak_rows, key=lambda x: np.nanmax(self.unnormalized_subpeak_rows)
            )
        return self.unnormalized_subpeak_rows
    #
    # def give_rows_of_secondary_subpeaks(self, coverage, chr_lens):
    #     if not hasattr(self, 'subpeak_rows'):
    #         self.give_rows_of_subpeaks_for_image(coverage, chr_lens)
    #     if len(self.subpeak_rows) > 1:
    #         return self.subpeak_rows[:-1]
    #     else:
    #         return []

    def give_row_of_highest_subpeak(self):
        if not hasattr(self, 'subpeak_rows'):
            self.give_rows_of_subpeaks_for_image(coverage, chr_lens)
        if len(self.subpeak_rows) > 0:
            return self.subpeak_rows[-1]
        else:
            return []

    def subpeak_stats(self, verbose=False):
        self.info = {'num_super_peaks': 1}
        self.info['num_peaks'] = len(self.subpeak_arr2d)
        assert isinstance(self.highest_subpeak, list)
        assert len(self.highest_subpeak) == 2  # (position, height)
        self.info['distances_to_highest_peak'] = [int(peak[0] - self.highest_subpeak[0]) for peak in self.subpeak_arr2d]
        self.info['absolute_distances_to_highest_peak'] = [abs(x) for x in self.info['distances_to_highest_peak']]
        stat_line = '''Super peak stats (iv: {iv})
        Number of peaks {np}
        Distances to highest peak {dp}
        Absolute distances to highest peak {ad}'''.format(
            iv=str(self.iv), np=self.info['num_peaks'],
            dp=self.info['distances_to_highest_peak'],
            ad=self.info['absolute_distances_to_highest_peak']
        )
        self.stats = _stats(self.info)
        if verbose: print stat_line
        return self.info

    def stats_of_a_2d_arr_of_subpeaks(self, unnormalized_subpeak_rows):
        max_vals = [(index, max(arr)) for index, arr in enumerate(unnormalized_subpeak_rows)]
        max_vals = sorted(max_vals, key=lambda x: x[1])
        highest = unnormalized_subpeak_rows[max_vals[-1][0]]
        secondary = []
        for index, val in max_vals:
            secondary.append(unnormalized_subpeak_rows[index])
        self.imageinfo = {'num_2d_arrays': 1}
        self.imageinfo['distances_to_highest'] = []

    def do_cwt(self, coverage_array):
        if max(coverage_array) < 2:
            return []
        peak_pos = []
        for a_pos in signal.find_peaks_cwt(
                coverage_array, np.arange(25, 325, 100), noise_perc=5):
            if (len(peak_pos) > 0) and (a_pos > (10 + peak_pos[-1])):
                peak_pos.append(a_pos)
            elif (len(peak_pos) == 0):
                peak_pos.append(a_pos)
        return peak_pos


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    return complement(s)


class _stats():

    def __init__(self, info, individual_dicts=None):
        assert type(info) is type({})
        self.info = info
        self.purge_zeroes(self.info)
        if individual_dicts is not None:
            self.individual_dicts = individual_dicts
        else:
            self.individual_dicts = [info]

    def purge_zeroes(self, adict):
        adict['distances_to_highest_peak'] = [val for val in adict['distances_to_highest_peak'] if val > 0]
        adict['absolute_distances_to_highest_peak'] = [val for val in adict['absolute_distances_to_highest_peak'] if val > 0]

    def __add__(self, y):
        #assert type(self) is type(y)
        self.purge_zeroes(y.info)
        self.individual_dicts.append(y.info)
        for key in self.info:
            if key in y.info:
                self.info[key] += y.info[key]
        #print 'after adding' + str(self.__dict__)
        return _stats(self.info, individual_dicts=self.individual_dicts)

    def __str__(self):
        self.calculate()
        return """
        num peaks {n}
        num super peaks {ns}
        num peaks per superpeak {pp}
        median peaks in superpeak {medp}
        median distance to highest peak {mh}
        median absolute dist to highest {mah}
        """.format(n=self.info['num_peaks'], ns=self.info['num_super_peaks'],
                   pp=self.peaks_per_superpeak, mh=self.median_dist_to_highest,
                   medp=self.median_peaks_in_superpeak,
                   mah=self.median_absolute_dist_to_highest,)

    def calculate(self):
        for x in [
                'distances_to_highest_peak',
                'absolute_distances_to_highest_peak']:
            self.info[x] = [float(y) for y in self.info[x]]
        for x in [
                'num_super_peaks',
                'num_peaks']:
            self.info[x] = float(self.info[x])
        self.median_peaks_in_superpeak = np.median(
            [int(x['num_peaks']) for x in self.individual_dicts])
        self.peaks_per_superpeak = (self.info['num_peaks']/self.info['num_super_peaks'])
        self.median_dist_to_highest = np.median(self.info['distances_to_highest_peak'])
        self.median_absolute_dist_to_highest = np.median(self.info['absolute_distances_to_highest_peak'])


class subpeak_rows():

    def __init__(self, coverage, from_array=None):
        if from_array is not None:
            self.coverage = from_array[0]
            self.rows = from_array
        else:
            assert (type(coverage) is type(np.array([]))) or (type(coverage) is type([]))
            if type(coverage) is type(np.array([])):
                coverage = coverage.tolist()
            self.rows = [coverage]
            self.coverage = coverage

    def __add__(self, y):
        self.rows.append(y.coverage)
        return subpeak_rows(self.coverage, from_array=self.rows)

    def __str__(self):
        if not hasattr(self, 'ave_y'):
            self.aves()
        return str(self.ave_y)

    def normalized_rows(self):
        self.norm_rows = []
        for row in sorted(self.rows, key=lambda x: max(x), reverse=True):
            try:
                maxval = float(max(row))
                self.norm_rows.append([float(i)/maxval for i in row])
            except:
                print "row:"
                print row
                sys.exit()
        return self.norm_rows

    def height_at_center(self, alist):
        return alist[int(len(alist)/2)]

    def sort_rows(self):
        self.rows = sorted(self.rows, key=lambda x: self.height_at_center(x))
        self.normalized_rows()
        return self.rows

    def aves(self):
        self.width = max([len(x) for x in self.rows])
        self.ave_y = [0 for x in range(0, self.width)]

        # Pad.
        for n, row in enumerate(self.rows):
            if len(row) < self.width:
                self.rows[n] = row + [0] * (self.width - len(row))
        for col in range(0, self.width):
            for row in self.rows:
                self.ave_y[col] += row[col]
            self.ave_y[col] = float(self.ave_y[col])/max([1., len(self.rows)])
        return self.ave_y

    def get_scaled_values(self):
        if not hasattr(self, 'ave_y'): self.aves()
        raw_integral = float(sum(self.ave_y))
        self.scaled_values = [y/raw_integral for y in self.ave_y]
        return self.scaled_values

    def fig(self, plt_object, line_params='k-'):
        if not hasattr(self, 'ave_y'): self.aves()
        self.get_scaled_values()
        plt_object.plot(range(0, len(self.ave_y)), self.scaled_values,  line_params)
        return plt_object


def add_subpeak_arr(subp_obj, to_add):
    if len(to_add) == 0:
        return subp_obj
    np = subpeak_rows(to_add)
    if subp_obj is None:
        subp_obj = np
    else:
        subp_obj = subp_obj + np
    return subp_obj


class sbGene:

    def __init__(self, iv):
        self.iv = iv

    def find_subpeaks(self, coverage, exons_dict):
        self.determine_coverage_accross_exons(coverage, exons_dict)
        self.peaks = self.do_cwt(self.exon_coverage)
        self.peak_ivs_genomic_coord = [HTSeq.GenomicPosition(
            self.iv.chrom,
            self.position_in_exon_coverage_mapped_to_genomic_position[pos],
            self.iv.strand,
            ) for pos in self.peaks]
        self.peak_ivs_relative = [HTSeq.GenomicPosition(
            self.iv.chrom, pos, self.iv.strand,
            ) for pos in self.peaks]

    def filter_by_intersection(self, list_of_ivs):
        self.filtered_peak_ivs_genomic_coord = []
        for iv in self.peak_ivs_genomic_coord:
            for alt_iv in list_of_ivs:
                if alt_iv[3] != iv.strand: continue
                if alt_iv[0] != iv.chrom: continue
                if alt_iv[1] <= iv.pos <= alt_iv[2]:
                    self.filtered_peak_ivs_genomic_coord.append(iv)
        return self.filtered_peak_ivs_genomic_coord

    def add_to_wig(self, ga):
        for iv in self.peak_ivs_genomic_coord:
            ga[iv] += 1

    def do_cwt(self, coverage_array):
        if max(coverage_array) < 2:
            return []
        peak_pos = []
        # for a_pos in signal.find_peaks_cwt(
        #         coverage_array, np.arange(15, 125, 20), noise_perc=5):
        #     if (len(peak_pos) > 0) and (a_pos > (1 + peak_pos[-1])) and (
        #             coverage_array[a_pos] >= 2):
        #         peak_pos.append(a_pos)
        #     elif (len(peak_pos) == 0):
        #         peak_pos.append(a_pos)
        try:
            arr_pos = find_local_extrema(coverage_array)[0]
        except:
            print "???"
            print coverage_array
            return []
        # print arr_pos
        for a_pos in arr_pos:
            if a_pos == 0: continue
            a = np.where(coverage_array==coverage_array[a_pos])
            height = coverage_array[a_pos]
            first_drop = np.argmax(coverage_array[a_pos:]<height) + a_pos
            right = np.max([first_drop - 1, a_pos])
            left = a_pos
            a_pos = int(float(left + right)/2.)
            # if np.max(coverage_array[max([1, left-10]):a_pos]) > coverage_array[a_pos]:
            #     continue
            # if np.max(coverage_array[a_pos+1:min([len(coverage_array), right+10])]) > coverage_array[a_pos]:
            #     continue
            if (len(peak_pos) > 0) and (a_pos > (10 + peak_pos[-1])) and (
                    coverage_array[a_pos] >= 2):
                peak_pos.append(a_pos)
            elif len(peak_pos) == 0:
                peak_pos.append(a_pos)
        return list(peak_pos[:])

    def determine_fbe_coverage_across_exons(self, sequences, exons_dict):
        self.exon_seq = ''
        for exon_num in sorted(exons_dict.keys(), key=lambda x: int(x)):
            (chrm, left, right, strand) = exons_dict[exon_num]
            iv = HTSeq.GenomicInterval(chrm, left, right, strand)
            self.exon_seq += get_seq_from_iv(iv, sequences)
        self.fbe_coverage = []
        last_pos = 0
        self.motif_no_star = '(' + self.motif + ')'
        for m in re.finditer(self.motif_no_star, self.exon_seq, re.IGNORECASE):
            dist_from_last = m.start() - last_pos
            if dist_from_last > 0:
                self.fbe_coverage.extend([0] * (dist_from_last))
            self.fbe_coverage.extend([1] * (m.end() - m.start()))
            last_pos = m.end()
        self.fbe_coverage.extend([0] * (len(self.exon_seq) - last_pos))
        self.fbe_coverage = np.array(self.fbe_coverage)
        return self.fbe_coverage

    def determine_coverage_accross_exons(
            self, coverage, exons_dict):
        '''
        exon    0      1       2
        genomic a      b       c
        rela    012345601234567012345
        '''
        self.exon_coverage = np.array([])
        self.position_in_exon_coverage_mapped_to_genomic_position = {}
        first_exon_left = exons_dict[sorted(exons_dict.keys(), key=lambda x: int(x))[0]][1]
        strand = '+'
        reverse = False
        if exons_dict.values()[0][-1] == '-':
            strand = '-'
            reverse = True
        for exon_num in sorted(exons_dict.keys(), key=lambda x: int(x), reverse=reverse):
            (chrm, left, right, strand) = exons_dict[exon_num]
            iv = HTSeq.GenomicInterval(chrm, left, right, strand)
            this_exon_coverage = np.fromiter(coverage[iv], dtype='f')
            for pos in range(
                    len(self.exon_coverage),
                    len(self.exon_coverage) + len(this_exon_coverage)):
                self.position_in_exon_coverage_mapped_to_genomic_position[pos] = left + pos - len(self.exon_coverage)
            self.exon_coverage = np.concatenate([
                self.exon_coverage, this_exon_coverage])
            sorted_pos = sorted(
                self.position_in_exon_coverage_mapped_to_genomic_position.keys(),
            key=lambda x: int(x))


def process_file(filename, coverage, txpts, sequences, chr_lens, args=None,
                 lib=None):
    print 'process_file():'
    output_dirname = 'figs/{d}/'.format(d=os.path.basename(os.path.dirname(filename)))
    print "Writing to {s}.".format(s=output_dirname)
    df = pandas.read_csv(filename, sep='\t')
    peaks = df.to_dict('records')
    gene_targets = list(set(df['gene_name'].tolist()))
    gene_targets = filter(lambda x: isinstance(x, str), gene_targets)
#    gene_targets = ['ced-12']
#    print txpts['ced-12']
#    sys.exit()
    suffix = ''
    if args is not None:
        suffix=args.use
    analyze_by_applying_cwt_to_all_reads_in_gene(
        gene_targets, coverage, txpts, sequences, chr_lens, peaks,
        output_dirname=output_dirname, lib=lib)
    analyze_by_applying_cwt_to_the_peak_regions(
        peaks, coverage, txpts, sequences, chr_lens, output_dirname=output_dirname,
        suffix=suffix, lib=lib)


def analyze_by_applying_cwt_to_all_reads_in_gene(
        gene_targets, coverage, txpts, sequences, chr_lens,
        called_peaks_list, output_dirname='figs/',
        lib=None):
    # Analysis by applying CWT to the total coverage in a gene.
    genes = []
    # subpeak_highest = None
    # subpeak_secondary = None
    import collections
    called_peaks_by_gene = collections.defaultdict(list)
    for peak_row in called_peaks_list:
        called_peaks_by_gene[peak_row['gene_name']].append(peak_row)
    subpeak_list = []
    subpeak_wig = HTSeq.GenomicArray('auto', stranded=True)
    for n, gene_name in enumerate(gene_targets):
        #if n > 100: break
        if gene_name not in txpts:
            print "Missing " + gene_name
            continue
        gene = txpts[gene_name]
#        if gene.strand == '-': continue
        sp = sbGene(HTSeq.GenomicInterval(
            gene.chrom, gene.txpt_left, gene.txpt_right, gene.strand))
        sp.gene_name = gene_name
        sp.motif = lib['motif']
#        sp.determine_coverage_accross_exons(coverage, txpts[gene_name].exons_dict)
        rc_corrected_exons_dict = {}
        if gene.strand == '-':
            num_exons = len(txpts[gene_name].exons_dict)
            for exon_num in txpts[gene_name].exons_dict:
                (chrm, left, right, strand) = txpts[gene_name].exons_dict[exon_num]
                _left = chr_lens[chrm] - right
                _right = chr_lens[chrm] - left
                rc_corrected_exons_dict[num_exons - exon_num + 1] = (chrm, _left, _right, '-')
                sp.find_subpeaks(coverage, rc_corrected_exons_dict)
        else:
            sp.find_subpeaks(coverage, txpts[gene_name].exons_dict)
        list_of_peaks_iv = []
        for peak in called_peaks_by_gene[gene_name]:
            list_of_peaks_iv.append(
                [peak['chrm'], peak['left'], peak['right'], peak['strand']])
        sp.peak_ivs_genomic_coord = sp.filter_by_intersection(list_of_peaks_iv)
        sp.add_to_wig(subpeak_wig)
        sp.determine_fbe_coverage_across_exons(sequences, txpts[gene_name].exons_dict)
        genes.append(sp)
        genomic_pos_list = []
        for peak_iv in sp.filtered_peak_ivs_genomic_coord:
            genomic_pos_list.append(
                "{chr}:{npos}:{strand}".format(
                    chr=peak_iv.chrom, npos=peak_iv.pos, strand=peak_iv.strand))
        subpeak_list.append({'peaks_per_gene': len(sp.filtered_peak_ivs_genomic_coord),
                             'gene_name': sp.gene_name,
                             'peak_positions': ", ".join(genomic_pos_list)})
    subpeak_wig.write_bedgraph_file('peaks_pos.wig', strand='+')
    subpeak_wig.write_bedgraph_file('peaks_neg.wig', strand='-')
    num_peaks_table = pandas.DataFrame(subpeak_list)
    print num_peaks_table
    num_peaks_table.sort(columns=['peaks_per_gene'], ascending=False, inplace=True)
    num_peaks_table.to_csv('subpeaks_per_gene_by_all_reads_in_gene.txt', sep='\t', index=False)
    num_peaks_table.to_csv(
        'subpeaks_per_gene.txt', sep='\t', index=False, columns=['gene_name', 'peaks_per_gene'])
    print "Converting to 2d list"
    list2d = list2d_of_coverage_across_txpts(genes)
    list2d_fbe = list2d_of_fbe_coverage_across_txpts(genes)
    image = build_image_objects_for_heatmaps.convert_to_image(list2d)
    print 'fbe:'
    image_fbe = build_image_objects_for_heatmaps.convert_to_image(
        list2d_fbe, normalize=False)
    print 'image fbe'
    print '*'
    for row in image:
        row[row < 0.01] = 0
    # Image is a 2d numpy array
    get_m_info(image)
    # image = build_image_objects_for_heatmaps.cluster_image(
    #     image,
    #     do_second_image=False,
    #     metric='correlation', false_color_part=None)
    outfname = output_dirname + '/across_txpts.pdf'
    print 'writing txpt heatmmaps to ' + outfname
    output_heatmap_figures.build_heatmap_from_arrays(
        None, image, output_filename=outfname,
        set_scale_from_data=True)
    output_heatmap_figures.build_heatmap_from_arrays(
        None, image_fbe, output_filename=outfname.partition('.pdf')[0] + '_fbes.pdf',
        set_scale_from_data=True)
    edge_of_txpt.make_graph_of_a_txpt_region(
        list2d, list2d_fbe, edge=200,
        outfname=output_dirname + '/coverage_in_region/coverage_in_region_of_txpt.pdf')

    # output_heatmap_figures.build_heatmap_from_one_array(
    #     image,
    #     false_color_part=None,
    #     output_filename=outfname,
    #     center_line=True)

    # Coverage of txpt is in sp.gene_arr
    # pulsar.make_fig(
    #     subpeak_highest.normalized_rows()[:70],
    #     output_filename=output_dirname + '/pulsar_highest_peak_in_gene.pdf')
    # pulsar.make_fig(
    #     subpeak_secondary.normalized_rows()[:70],
    #     output_filename=output_dirname + '/pulsar_secondary_peaks_in_gene.pdf')
    # subpeak_highest.fig(output_dirname + '/average_secondary_peak_in_gene.pdf')
    # subpeak_secondary.fig(output_dirname + '/average_highest_peak_in_gene.pdf')


def list2d_of_coverage_across_txpts(genes):
    assert isinstance(genes, list)
    if len(genes) < 1:
        return None
    assert isinstance(genes[0], sbGene)
    list2d = []
    for gene in genes:
        list2d.append(
            normalize_row_to_highest(gene.exon_coverage.tolist(), as_list=True)
        )
    list2d = sorted(list2d, key=lambda x: len(x))
    return list2d


def list2d_of_fbe_coverage_across_txpts(genes):
    assert isinstance(genes, list)
    if len(genes) < 1:
        return None
    assert isinstance(genes[0], sbGene)
    list2d = []
    for gene in genes:
        list2d.append(
            gene.fbe_coverage.tolist()
        )
    list2d = sorted(list2d, key=lambda x: len(x))
    return list2d


def get_seq_from_iv(_iv, sequences):
    if _iv.end > len(sequences[_iv.chrom]):
        end = len(sequences[_iv.chrom])
    else: end = _iv.end
    start = max([_iv.start, 0])
    seq = sequences[_iv.chrom][start:end]
    if _iv.strand == '-':
        seq = rc(seq)
    return seq


def analyze_by_applying_cwt_to_the_peak_regions(
        peaks, coverage, txpts, sequences, chr_lens, output_dirname='figs/', suffix='',
        lib=None):
    print 'Analyzing by applying cwt to unmerge peaks in peak regions.'
    # Analysis by superpeak.
    supers = []
    for row in peaks:
        sp = superpeak(HTSeq.GenomicInterval(
            row['chrm'], row['left'], row['right'], row['strand']))
        sp.motif = lib['motif']
        sp.find_subpeaks(coverage)
        if sp.highest_subpeak is None or len(sp.subpeak_arr2d) == 0:
            print "No peaks found in %s..." % row['gene_name']
            continue
        sp.gene_name = row['gene_name']
        supers.append(sp)
    for peak in supers:
        peak.subpeak_stats()
        peak.get_seqs_from_ivs(sequences)
    write_fastas(supers, sequences, suffix=suffix, motif=lib['motif'])
    score_fastas(suffix=suffix, motif=lib['motif'])
    how_many_genes_have_secondary_peaks(supers)
    pulsar_of_a_subset_of_peaks(supers, coverage, chr_lens,
                                suffix=suffix, output_dirname=output_dirname)
    stats_on_all_superpeaks(supers)
    make_pulsars_and_average_peak_graphs(
        supers, coverage, chr_lens, output_dirname=output_dirname, suffix=suffix)
    make_image_array_and_draw_heatmaps(
        supers, coverage, chr_lens, output_dirname=output_dirname, suffix=suffix,
                     output_basename='subpeaks_heatmap.pdf')


def pulsar_of_a_subset_of_peaks(superpeak_list, coverage, chr_lens,
                                output_dirname='figs/', suffix='',
                                names=None, perspective=True,
                                two_color=None):
    genes_with_high_subpeak_numbers = [
        'plk-3', 'gck-1', 'math-33', 'sel-8', 'gld-2', 'mpk-1', 'ZK616.5', 'ifet-1',
        'F53C11.4', 'T04C4.1']
    superpeaks_by_gene = list_to_dict(superpeak_list)
    subset_superpeak_list = []
    by_num = sorted(superpeaks_by_gene.keys(), key=lambda x: len(superpeaks_by_gene[x]))
    # for gene in genes_with_high_subpeak_numbers:
    #     if gene in superpeaks_by_gene:
    #         subset_superpeak_list.extend(superpeaks_by_gene[gene])
    #     else:
    #         print "Not a target: %s" % gene
    for x in by_num[-10:]:
        subset_superpeak_list.extend(superpeaks_by_gene[x])
    suffix += 'most_peaks'
    make_pulsars_and_average_peak_graphs(
        subset_superpeak_list, coverage, chr_lens, output_dirname=output_dirname,
        suffix=suffix,  names=names, perspective=perspective,
    two_color=two_color)


def how_many_genes_have_secondary_peaks(superpeak_list):
    superpeaks_by_gene = list_to_dict(superpeak_list)
    num_peaks_per_gene = {}
    abs_distances_to_primary_peak = []
    distances_to_primary_peak = []
    ninc = 0
    num_secondary_by_pos_list = 0
    for gene in superpeaks_by_gene:
        pos_and_height_of_each_subpeak = []
        num_peaks_per_gene[gene] = 0
        highest_subpeak = (0, 0)
        for peak in superpeaks_by_gene[gene]:
            if peak.highest_subpeak[1] > highest_subpeak[1]:
                highest_subpeak = peak.highest_subpeak
        for peak in superpeaks_by_gene[gene]:
            num_peaks_per_gene[gene] += len(peak.subpeak_arr2d)
            for subpeak in peak.subpeak_arr2d:
                pos_and_height_of_each_subpeak.append([subpeak[0], subpeak[1]])
        num_secondary_by_pos_list += max([
            len(pos_and_height_of_each_subpeak) - 1, 0])
        if len(pos_and_height_of_each_subpeak) > 1:
            for x in pos_and_height_of_each_subpeak[:-1]:
                if (x[0] - highest_subpeak[0]) == 0:
                    ninc += 1
                    #print 'uh oh'
            abs_distances_to_primary_peak.extend([
                abs(x[0] - highest_subpeak[0]) for x in pos_and_height_of_each_subpeak[:-1] if x[0] != highest_subpeak[0]])
            distances_to_primary_peak.extend([
                int(x[0] - highest_subpeak[0]) for x in pos_and_height_of_each_subpeak[:-1] if x[0] != highest_subpeak[0]])
    print 'ninc'
    print ninc
    print """
    Intrapeak distances (to primary peak) measured {np}
    Number of secondary peaks observed {nsl}
    Number of genes {ng}
    Median absolute distance to primary peak {medd}
    Mean absolute distance to primary peak {meand}
    Median relative distance to primary peak {redmedd}
    Mean relative distance to primary peak {redmean}
    """.format(
        np=len(distances_to_primary_peak),
        nsl=num_secondary_by_pos_list,
        ng=len(superpeaks_by_gene),
        medd=np.median(abs_distances_to_primary_peak),
        meand=np.mean(abs_distances_to_primary_peak),
        redmedd=np.median(distances_to_primary_peak),
        redmean=np.mean(distances_to_primary_peak),
    )
    by_peak_number = sorted(
        num_peaks_per_gene.keys(),
        key=lambda x: num_peaks_per_gene[x])
    if len(by_peak_number) > 9:
        print "Genes with the most peaks:"
        for n in range(1, 11):
            print "{rank}.: {gene}\t{peak_number}".format(
                rank=n, gene=by_peak_number[-n],
                peak_number=num_peaks_per_gene[by_peak_number[-n]]
            )
        print "Genes with the fewest peaks:"
        for n in range(0, 10):
            print "{rank}.: {gene}\t{peak_number}".format(
                rank=n, gene=by_peak_number[n],
                peak_number=num_peaks_per_gene[by_peak_number[n]]
            )
    peaks_per_gene_table = pandas.DataFrame(
        zip(num_peaks_per_gene.keys(), num_peaks_per_gene.values()),
        columns=['gene_name', 'peaks_per_gene']
    )
    peaks_per_gene_table.sort(columns='peaks_per_gene', inplace=True,  ascending=False)
    peaks_per_gene_table.to_csv('subpeaks_per_gene_by_only_reads_in_peak_regions.txt',
                                sep='\t', index=False)
    num_with_secondary_peaks = len([x for x in num_peaks_per_gene.values() if x>1])
    num_without_secondary_peaks = len([x for x in num_peaks_per_gene.values() if x==1])
    num_without_peaks = len([x for x in num_peaks_per_gene.values() if x<1])
    print """
    Gene targets {tg}
    Genes with secondary peaks {nw}
    Genes without secondary peaks {nwithout} ({perc}%)
    Genes without peaks {nop}""".format(
        tg=len(superpeaks_by_gene),
        nw=num_with_secondary_peaks,
        nwithout=num_without_secondary_peaks,
        perc="%.3f" % float(100. * float(num_without_secondary_peaks)/float(len(superpeaks_by_gene))),
        nop=num_without_peaks,
    )


def score_fastas(suffix='', motif='.*[tT][gG][tT]\w\w\w[Aa][Tt].*'):
    with open('data/fasta/subpeaks/highest_%s.fa' % suffix, 'r') as f:
        seqs_h = []
        for li in f:
            if re.match('\A>.*', li) is None:
                seqs_h.append(li.rstrip('\n'))
    frac_in_primary = fraction_of_seqs_with_fbe(seqs_h, motif=motif)
    with open('data/fasta/subpeaks/secondary_%s.fa' % suffix, 'r') as f:
        seqs_l = []
        for li in f:
            if re.match('\A>.*', li) is None:
                seqs_l.append(li.rstrip('\n'))
    frac_in_secondary = fraction_of_seqs_with_fbe(seqs_l, motif=motif)
    print '''
    score_fastas():
    Primary peaks {np}. Secondary peaks {ns}.
    Fraction of primary peaks with FBE {fp}.
    Fraction of secondary peaks with FBE {fs}.
    Median length of primary peak {medpl}
    Median length of secondary peak {secpl}'''.format(
        np=len(seqs_h), ns=len(seqs_l),
        fp="%.3f" % float(100. * frac_in_primary),
        fs="%.3f" % float(100. * frac_in_secondary),
        medpl=np.median([len(x) for x in seqs_h]),
        secpl=np.median([len(x) for x in seqs_l]),
    )


def fraction_of_seqs_with_fbe(seq_list, motif=None):
    with_fbe = 0
    for seq in seq_list:
        if re.match(motif, seq, re.IGNORECASE) is not None:
            with_fbe += 1
    return float(with_fbe)/float(max(1, len(seq_list)))


def get_iv_of_peaks_sorted_from_superpeaks_list(
        superpeak_list, width=35, get_txpt_bounds=None):
    superpeaks_by_gene = list_to_dict(superpeak_list)
    peak_iv_by_gene = {}  # Dict of sorted lists, highest peak last.
    for gene in superpeaks_by_gene:
        peak_iv_by_gene[gene] = []
        pos_and_height_of_each_subpeak = []
        highest_subpeak = (0, 0)
        for peak in superpeaks_by_gene[gene]:
            if peak.highest_subpeak[1] > highest_subpeak[1]:
                highest_subpeak = peak.highest_subpeak
        for peak_index, peak in enumerate(superpeaks_by_gene[gene]):
            for subpeak in peak.subpeak_arr2d:
                pos_and_height_of_each_subpeak.append([subpeak[0], subpeak[1], peak_index])
        pos_and_height_of_each_subpeak = sorted(
            pos_and_height_of_each_subpeak,
            key=lambda x: x[1]
        )
        for pos, height, peak_index in pos_and_height_of_each_subpeak:
            left = superpeaks_by_gene[gene][peak_index].iv.start
            iv = HTSeq.GenomicInterval(
                superpeaks_by_gene[gene][peak_index].iv.chrom,
                max([left + pos - width, 1]), left + pos + width,
                superpeaks_by_gene[gene][peak_index].iv.strand,
            )
            peak_iv_by_gene[gene].append(iv)
    return peak_iv_by_gene


def write_fastas(superpeak_list, sequences, suffix='',
                 motif='.*[tT][gG][tT]\w\w\w[aA][tT].*'):
    superpeaks_by_gene = list_to_dict(superpeak_list)
    seqs = []
    seqs_under_highest = []
    seqs_under_secondary = []
    peak_iv_by_gene = get_iv_of_peaks_sorted_from_superpeaks_list(
        superpeak_list, width=35)
    for gene in peak_iv_by_gene:
        gene_seqs = []
        iv_list = peak_iv_by_gene[gene]
        for iv in iv_list:
            gene_seqs.append(get_seq_from_iv(iv, sequences))
        seqs_under_highest.append(gene_seqs[-1])
        seqs_under_secondary.extend(gene_seqs[:-1])
        seqs.extend(gene_seqs)
    high = ''
    high_wo_fbe = ''
    high_w_fbe = ''
    for n, seq in enumerate(seqs_under_highest):
        high += '>{n}\n{s}\n'.format(n=n, s=seq)
        if re.match(motif, seq, re.IGNORECASE) is not None:
            high_w_fbe += '>{n}\n{s}\n'.format(n=n, s=seq)
        else:
            high_wo_fbe += '>{n}\n{s}\n'.format(n=n, s=seq)
    low = ''
    low_wo_fbe = ''
    low_w_fbe = ''
    for n, seq in enumerate(seqs_under_secondary):
        low += '>{n}\n{s}\n'.format(n=n, s=seq)
        if re.match(motif, seq, re.IGNORECASE) is not None:
            low_w_fbe += '>{n}\n{s}\n'.format(n=n, s=seq)
        else:
            low_wo_fbe += '>{n}\n{s}\n'.format(n=n, s=seq)
    if not os.path.exists('data/fasta/subpeaks'):
        os.system('mkdir data/fasta/subpeaks')
    with open('data/fasta/subpeaks/highest_%s.fa' % suffix, 'w') as f:
        f.write(high)
    with open('data/fasta/subpeaks/secondary_%s.fa' % suffix, 'w') as f:
        f.write(low)
    with open('data/fasta/subpeaks/highest_w_fbe_%s.fa' % suffix, 'w') as f:
        f.write(high_w_fbe)
    with open('data/fasta/subpeaks/secondary_w_fbe_%s.fa' % suffix, 'w') as f:
        f.write(low_w_fbe)
    with open('data/fasta/subpeaks/highest_wo_fbe_%s.fa' % suffix, 'w') as f:
        f.write(high_wo_fbe)
    with open('data/fasta/subpeaks/secondary_wo_fbe_%s.fa' % suffix, 'w') as f:
        f.write(low_wo_fbe)


def stats_on_all_superpeaks(superpeaks):
    assert type(superpeaks) is type([])
    # Simple stats.
    all_subpeak_stats = superpeaks[0].stats
    for peak in superpeaks[1:]:
        all_subpeak_stats = all_subpeak_stats + peak.stats
    print "By superpeak:"
    print str(all_subpeak_stats)


def image_from_iv_list(iv_list, coverage):
    subpeak_rows = []
    for subpeak_iv in iv_list:
        unnormalized_subpeak_row = np.nan_to_num(np.fromiter(coverage[subpeak_iv], dtype='i'))
        unnormalized_subpeak_row = np.clip(unnormalized_subpeak_row, 1e-9, 1e9)
        subpeak_rows.append(normalize_row_to_highest(unnormalized_subpeak_row, as_list=True))
    return subpeak_rows


def make_pulsars_and_average_peak_graphs(
        superpeak_list, coverage, chr_lens, output_dirname='figs/',
        suffix='', names=None, perspective=None,
        two_color=None):
    """
    A list of superpeaks is passed, which is reorganized to a dict by gene.
    The callpeaks identified peaks are investigated.
    """
    # Imaging the averages.
    # Collect data into a by-gene 2d array and a genome 2d array.
    # Specifically, we build arrays of signal around all subpeaks.
    peak_iv_by_gene = get_iv_of_peaks_sorted_from_superpeaks_list(
        superpeak_list, width=1000)
    highest_peak_ivs = []
    highest_peak_names = []
    secondary_peak_ivs = []
    secondary_peak_names = []
    for gene in peak_iv_by_gene:
        if len(peak_iv_by_gene) > 0:
            highest_peak_ivs.append(peak_iv_by_gene[gene][-1])
            highest_peak_names.append(gene)
            secondary_peak_ivs.extend(peak_iv_by_gene[gene][:-1])
            secondary_peak_names.extend([gene] * len(peak_iv_by_gene[gene][:-1]))
    highest_image = image_from_iv_list(highest_peak_ivs, coverage)
    secondary_image = image_from_iv_list(secondary_peak_ivs, coverage)
    # Pulsar accepts lists or np arrays.
    num_rows_highest = min([70, len(highest_image)])
    num_rows_secondary = min([70, len(secondary_image)])
    pulsar.make_fig(
        highest_image[:num_rows_highest],
        output_filename=output_dirname + '/pulsar_highest_peak_%s.pdf' % suffix,
        names=highest_peak_names, perspective=perspective,
        two_color=two_color)
    pulsar.make_fig(
        secondary_image[:num_rows_secondary],
        output_filename=output_dirname + '/pulsar_secondary_peaks_%s.pdf' % suffix,
        names=secondary_peak_names, perspective=perspective,
        two_color=two_color)
    highest_subpeak_in_gene_rows = subpeak_rows(None, from_array=highest_image)
    secondary_subpeak_rows = subpeak_rows(None, from_array=secondary_image)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    plt.clf()
    plt.xlabel('Nucleotides')
    plt.ylabel('Average coverage in a peak')
    plt.rcParams.update({'font.size': 8})
    fig.set_size_inches(4, 4)
    plt = highest_subpeak_in_gene_rows.fig(plt, line_params='k-')
    plt = secondary_subpeak_rows.fig(plt, line_params='r:')
    max_y = max([
        max(highest_subpeak_in_gene_rows.get_scaled_values()),
        max(secondary_subpeak_rows.get_scaled_values())])
    graph_y_max = 1.1 * max_y
    plt.plot([950, 950], [0, graph_y_max], 'k--')
    plt.plot([1050, 1050], [0, graph_y_max], 'k--')
    plt.xticks([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800],
               ['-800', '-600', '-400', '-200', '0',
                '200', '400', '600', '800'])
        #[0, 200, 400, 600, 800, 1e3])
    plt.ylim([0, graph_y_max])
    plt.xlim([600, 1400])
    fname = output_dirname + '/average_highest_peak_%s.pdf' % suffix
    plt.savefig(fname, format='pdf')


def list_to_dict(superpeaks):
    import collections
    assert type(superpeaks) is type([])
    by_gene = collections.defaultdict(list)
    for sp in superpeaks:
        # by_gene.setdefault(sp.gene_name, [])
        by_gene[sp.gene_name].append(sp)
    return by_gene


def make_image_array_and_draw_heatmaps(superpeak_list, coverage, chr_lens,
                     output_dirname='figs/methods/', output_basename='subpeaks_heatmap.pdf',
                                       suffix=''):
    """
    """
    #os.makedirs(output_dirname)
#    superpeaks_by_gene = list_to_dict(superpeaks)
    # Collect data into a by-gene 2d array and a genome 2d array.
    # Specifically, we build arrays of signal around all subpeaks.
    peak_iv_by_gene = get_iv_of_peaks_sorted_from_superpeaks_list(
        superpeak_list, width=1e3)
    highest_peak_ivs = []
    all_peak_ivs = []
    secondary_peak_ivs = []
    for gene in peak_iv_by_gene:
        if len(peak_iv_by_gene) > 0:
            highest_peak_ivs.append(peak_iv_by_gene[gene][-1])
            secondary_peak_ivs.extend(peak_iv_by_gene[gene][:-1])
            all_peak_ivs.extend(peak_iv_by_gene[gene][:])
    list2d_highest_peak = image_from_iv_list(highest_peak_ivs, coverage)
    list2d_secondary_peaks = image_from_iv_list(secondary_peak_ivs, coverage)
    list2d_all_peaks = image_from_iv_list(all_peak_ivs, coverage)
    convert_to_image_and_write(
        list2d_all_peaks,
        output_filename=output_dirname + '/subpeaks_all_peaks_%s.pdf' % suffix)
    convert_to_image_and_write(
        list2d_highest_peak,
        output_filename=output_dirname + '/subpeaks_highest_peaks_%s.pdf' % suffix)
    convert_to_image_and_write(
        list2d_secondary_peaks,
        output_filename=output_dirname + '/subpeaks_secondary_peaks_%s.pdf' % suffix)


def convert_to_image_and_write(list2d, output_filename='peaks_heatmap.pdf'):
#    print 'writing to {d}'.format(d=output_filename)
    image = build_image_objects_for_heatmaps.convert_to_image(list2d)
    # Image is a 2d numpy array
    get_m_info(image)
    image = build_image_objects_for_heatmaps.cluster_image(
        image,
        do_second_image=False,
        metric='euclidean', false_color_part=None)
    output_heatmap_figures.build_heatmap_from_one_array(
        image,
        false_color_part=None,
        output_filename=output_filename,
        center_line=True)


def get_m_info(m):
    size_range = set()
    types = set()
    val_types = set()
    for p in m:
        size_range |= set([len(p)])
        types |= set([type(p)])
        for v in p:
            val_types |= set([type(v)])


def normalize_row_to_highest(row, as_list=False):
    highest_val = float(np.max(row))
    if highest_val <= 0.:
        return row
    row = np.array([float(x)/highest_val for x in row])
    if as_list:
        return row.tolist()
    return row


def get_bedgraph(
        do_combine_bedgraphs=False, bedgraphs_folder='data/wigs_five_prime/',
        lib=None, args=None):
    if args is not None and args.use == 'lib':
        plus_file = lib['bedgraph_exp_plus']
        minus_file = lib['bedgraph_exp_minus']
        ga = {}
        ga['both'] = get_a_bedgraph(plus_file, minus_file)
        return ga
    plus_file = lib['bedgraph_exp_plus']
    minus_file = lib['bedgraph_exp_minus']
    ga = {}
    ga['both'] = get_a_bedgraph(plus_file, minus_file)
    if args is not None:
        if args.use == 'fbf1':
            plus_file = lib['bedgraph_fbf1_plus']
            minus_file = lib['bedgraph_fbf1_minus']
            ga['fbf1'] = get_a_bedgraph(plus_file, minus_file)
        elif args.use == 'fbf2':
            plus_file = lib['bedgraph_fbf2_plus']
            minus_file = lib['bedgraph_fbf2_plus']
            ga['fbf2'] = get_a_bedgraph(plus_file, minus_file)
    return ga


def get_a_bedgraph(plus_file, minus_file):
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    with open(plus_file, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] = float(s[3])
    with open(minus_file, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] = float(s[3])
    return ga


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s


def get_sequences(lib):
    fasta_filename = '/scratch/indexes/WS235.fa'
    sequences = dict(
        (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
    rc_sequences = dict(
        (p.name.split(' ')[0], rc(p.seq)) for p in HTSeq.FastaReader(fasta_filename))
    #chr_lens = dict(
    #    [(name, len(sequences[name])) for name in sequences])
    chr_lens = {}
    with open (lib['chr_sizes'], 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            chr_lens[s[0]] = int(s[1])
    return (sequences, rc_sequences, chr_lens)


def run(
        peak_fname, ga, lib, txpts, sequences, rc_sequences,
        chr_lens, args=None):
    if peak_fname is None:
        peak_fname = 'methods/filtered_gauss_nil_01_five_reps/combined_fbf.txt'
    logger.info('%s: Reading sequence information.' % datetime.datetime.now().strftime('%Hh%Mm'))
    if args is None:
        process_file(peak_fname, ga['both'], txpts, sequences, chr_lens, lib=lib)
    else:
        if args.use == 'fbf1':
            print "Processing fbf1 only."
            process_file(peak_fname, ga['fbf1'], txpts, sequences, chr_lens, args=args,
                         lib=lib)
        elif args.use == 'fbf2':
            print 'Processing fbf2 only.'
            process_file(peak_fname, ga['fbf2'], txpts, sequences, chr_lens, args=args,
                         lib=lib)
        else:
            process_file(peak_fname, ga['both'], txpts, sequences, chr_lens, args=args,
                         lib=lib)


if __name__ == '__main__':
    peak_fname = sys.argv[1]
    ga = get_bedgraph(
        bedgraphs_folder='/groups/Kimble/Aman Prasad/clip/data/wigs_coverage/')
    process_file(peak_fname, ga['both'])
