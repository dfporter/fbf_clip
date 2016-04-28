__author__ = 'dp'
import HTSeq
import numpy as np


class sbGene:

    def __init__(self, iv):
        self.iv = iv

    def find_subpeaks(self, coverage):
        width = 0
        expanded_iv = HTSeq.GenomicInterval(
            self.iv.chrom,
            max([self.iv.start - width, 1]),
            self.iv.end + width, self.iv.strand,
        )
        self.gene_arr = np.fromiter(coverage[expanded_iv], dtype='i')
        self.peaks = self.do_cwt(self.peak_arr)
        self.peak_ivs_relative = [HTSeq.GenomicPosition(
            self.iv.chrom, pos, self.iv.strand,
            ) for pos in self.sub_peaks]
        self.subpeak_arr2d = [[pos, self.gene_arr[pos]] for pos in self.peaks]
        self.subpeak_arr2d = sorted(self.subpeak_arr2d, key=lambda x: x[1])
        self.highest_subpeak = self.subpeak_arr2d[-1]
        self.subpeak_stats()

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

    def subpeak_stats(self, verbose=False):
#        return stats_of_a_2d_arr_of_subpeaks(self.subpeak_arr2d)
        self.info = {'num_super_peaks': 1}
        self.info['num_peaks'] = len(self.subpeak_arr2d)
        assert type(self.highest_subpeak) is type([])
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
