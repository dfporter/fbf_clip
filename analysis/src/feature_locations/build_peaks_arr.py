"""
Produces dict objects by (key=arbitrary number) of numpy arrays of peaks.

Takes as input:
    peaks: a pandas DataFrame
    txpts: a dict of flocs objects

Outputs:
    locsarr: dict of numpy arrays holding the CDS/txpt ranges for a gene
    peakslocs: dict of numpy arrays holding peak locations for a gene
    utrsarr: dict of numpy arrays holding the 3'UTR ranges for a gene
    peaks_in_utrs: dict of numpy arrays holding the peak locations in 3'UTRs

"""
__author__ = 'dp'
import pandas
import numpy as np
import logging
import datetime
logger = logging.getLogger(__name__)


def build_intrapeak_distances(peaks, txpts):
    logger.info('%s: Build intrapeak distances called.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    logger.info('peaks {n}. transcripts {p}'.format(
        n=len(peaks.index), p=len(txpts)
    ))
    locs = []
    for g in txpts:
        locs.append(len(txpts[g].peak_locs))
    logger.info("Average number of peak locs {mn}".format(mn=np.mean(locs)))
    n = -1
    all_intra_peak_distances = []
    intra_peak_distances_by_gene = {}
    genes_with_peaks = set()
    genes_with_multi = set()
    for gene_name in set(peaks['gene_name']):
        n += 1
        if gene_name not in txpts: continue
        if len(txpts[gene_name].peak_locs) == 0: continue
        genes_with_peaks |= set(gene_name)
        txpts[gene_name].get_distances_between_peaks()
        if len(txpts[gene_name].intrapeak_distances):
            genes_with_multi |= set(gene_name)
        for dist in txpts[gene_name].intrapeak_distances:
            all_intra_peak_distances.append(dist)
    logger.info('Genes with peaks {n}. Genes with peak pairs {p}'.format(
        n=len(list(genes_with_peaks)), p=len(list(genes_with_multi))
    ))
    return all_intra_peak_distances


def build_peak_raw_utr_signal(peaks, txpts, include_motif=False):
    utrs_dict = {}
    peaks_in_utrs_dict = {}
    fbes_in_utrs_dict = {}
    fbes_in_txpt_dict = {}
    n = -1
    missing = []
    for gene_name in set(peaks['gene_name']):
        n += 1
        if gene_name not in txpts:
            missing.append(gene_name)
            continue
        if include_motif:
            _fbe_in_txpt = [1] * (10 + len(txpts[gene_name].seq))
            for index, a_fbe in enumerate(txpts[gene_name].motif_locs):
                for pos in range(a_fbe[0], a_fbe[1]):
                    _fbe_in_txpt[pos] = 2
        if gene_name not in txpts: continue
        hits = peaks[peaks['gene_name']==gene_name]
        if include_motif:
            fbes_in_txpt_dict[n] = _fbe_in_txpt
        #if """CDS""" not in hits['location'].tolist():
        #    continue
        if len(txpts[gene_name].peak_locs) == 0: continue
        start_pos = int(txpts[gene_name].start_pos_in_seq)
        stop_pos = int(txpts[gene_name].stop_pos_in_seq)
        if (stop_pos - start_pos) < (0.3 * len(txpts[gene_name].seq)):
            continue
        _locsarr = [0.3] * len(txpts[gene_name].seq)
        utr_len =len(
            range(
                min(len(_locsarr), stop_pos + 1), len(_locsarr), 1)
        )
        _utr = [0.3] * len(txpts[gene_name].utr_arr)
        #_peak = [np.nan] * len(txpts[gene_name].utr_arr)
        if include_motif:
            _fbe = [1] * len(txpts[gene_name].utr_arr)
            for index, a_fbe in enumerate(txpts[gene_name].motif_locs):
                if a_fbe[0] <= txpts[gene_name].stop_pos_in_seq:
                    continue
                left = a_fbe[0] - txpts[gene_name].stop_pos_in_seq
                right = a_fbe[1] - txpts[gene_name].stop_pos_in_seq
                if right > len(txpts[gene_name].utr_arr):
                    continue
                for pos in range(left, right):
                    _fbe[pos] = 2
        if len(txpts[gene_name].utr_arr) < 5:
            continue
        utrs_dict[n] = _utr
        if include_motif:
            fbes_in_utrs_dict[n] = _fbe
        peaks_in_utrs_dict[n] = txpts[gene_name].utr_arr
    logger.warn("Missing in txpts (<- gtf) object: %s." % ", ".join(
        [str(x) for x in missing]))
    if include_motif:
        return (utrs_dict, peaks_in_utrs_dict, fbes_in_utrs_dict)
    else:
        return (utrs_dict, peaks_in_utrs_dict) 


def build_peak_raw_peak_region(peaks, txpts):
    highest_peak_raw = []
    secondary_peaks_raw = []
    secondary_peaks_highest_marked = []
    region_around_highest_raw = []
    highest_if_have_secondary_raw = []
    highest_if_have_secondary_highest_marked = []
    all_peaks_raw = []
    all_peaks_highest_marked = []
    region_around_all_peaks_raw = []
    for n, gene_name in enumerate(set(peaks['gene_name'])):
        if gene_name not in txpts: continue
        if len(txpts[gene_name].peak_locs) == 0: continue
        for _peak in txpts[gene_name].all_peaks_raw:
            all_peaks_raw.append(np.array(_peak))
        for _marks in txpts[gene_name].all_peaks_highest_marked:
            all_peaks_highest_marked.append(np.array(_marks))
        for _peak in txpts[gene_name].region_around_all_peaks_raw:
            region_around_all_peaks_raw.append(
                [np.array(_peak[0]), np.array(_peak[1]),
                 int(np.nanmax(_peak[0] + _peak[1]))]
            )
        highest_peak_raw.append(txpts[gene_name].highest_peak_raw)
        if len(txpts[gene_name].secondary_peaks_raw) > 0:
            highest_if_have_secondary_raw.append(txpts[gene_name].highest_peak_raw)
        if len(txpts[gene_name].secondary_peaks_raw) > 0:
            for _mark in txpts[gene_name].highest_peak_highest_marked:
                highest_if_have_secondary_highest_marked.append(np.array(_mark))
        for _peak in txpts[gene_name].secondary_peaks_raw:
            secondary_peaks_raw.append(np.array(_peak))
        for _mark in txpts[gene_name].secondary_peaks_highest_marked:
            secondary_peaks_highest_marked.append(np.array(_mark))
        _region = txpts[gene_name].region_around_highest_raw
        region_around_highest_raw.append(
            [np.array(_region[0]), np.array(_region[1]),
             int(np.nanmax(_region[0] + _region[1]))]
        )
    dont_sort = '''
    highest_peak_raw = sorted(highest_peak_raw, key=lambda x: np.nanmax(x))
    highest_if_have_secondary_raw = sorted(
        highest_if_have_secondary_raw, key=lambda x: np.nanmax(x))
    highest_if_have_secondary_highest_marked = np.array(
        sorted(highest_if_have_secondary_highest_marked, key=lambda x: np.nanmax(x))
    )
    '''
    highest_if_have_secondary_highest_marked = np.array(
        highest_if_have_secondary_highest_marked
    )
    all_peaks_highest_marked = np.array(all_peaks_highest_marked)
    secondary_peaks_highest_marked = np.array(secondary_peaks_highest_marked)
    dont_sort = '''
    highest_peak_raw = sorted(highest_peak_raw, key=lambda x: np.nanmax(x))
    all_peaks_raw = sorted(all_peaks_raw, key=lambda x: np.nanmax(x))
    all_peaks_highest_marked = np.array([
        x for (y, x) in sorted(
            zip(all_peaks_raw, all_peaks_highest_marked),
            key=lambda tup: np.nanmax(tup[0]))
    ])
    secondary_peaks_raw = sorted(secondary_peaks_raw, key=lambda x: np.nanmax(x))
    secondary_peaks_highest_marked = [
        x for (y, x) in sorted(
            zip(secondary_peaks_raw, secondary_peaks_highest_marked),
            key=lambda tup: np.nanmax(tup[0]))
    ]
    region_around_all_peaks_raw = sort_regions_and_reshape(region_around_all_peaks_raw)
    region_around_highest_raw = sort_regions_and_reshape(region_around_highest_raw)
    '''
    logger.info(
        '\t\tNumber of columns of all_peaks_raw (not as array) %i' % len(highest_peak_raw))
    return (all_peaks_raw, region_around_all_peaks_raw, highest_peak_raw,
            secondary_peaks_raw, region_around_highest_raw,
            all_peaks_highest_marked,
            secondary_peaks_highest_marked, highest_if_have_secondary_raw,
            highest_if_have_secondary_highest_marked)


def sort_regions_and_reshape(regions):
    regions = sorted(regions, key=lambda x: x[2])
    return [
        [x[0] for x in regions],
        [x[1] for x in regions]]


def build_peak_dicts(peaks, txpts):
    '''Create peak, txpt and utr signal objects.
    Returns:
    '''
    txpts_dict = {}
    peaks_dict = {}
    utrs_dict = {}
    peaks_in_utrs_dict = {}
    n = -1
    for gene_name in set(peaks['gene_name']):
        n += 1
        if gene_name not in txpts: continue
        if len(txpts[gene_name].peak_locs) == 0: continue
        start_pos = int(txpts[gene_name].start_pos_in_seq)
        stop_pos = int(txpts[gene_name].stop_pos_in_seq)
        if (stop_pos - start_pos) < (0.3 * len(txpts[gene_name].seq)):
            continue
        _locsarr = [0.3] * len(txpts[gene_name].seq)
        _peakslocs = [0] * len(txpts[gene_name].seq)
        for i in range(max(0, start_pos), min(len(_locsarr), stop_pos + 1)):
            _locsarr[i] = 0.7
        utr_len =len(
            range(
                min(len(_locsarr), stop_pos + 1), len(_locsarr), 1)
        )
        _utrsarr = [1.0] * utr_len
        _peaks_in_utr = [0] * utr_len
        if has_peak_in_range(txpts, gene_name, len(_peakslocs)):
            txpts_dict[n] = _locsarr
            add_peak_signal(txpts, gene_name, _peakslocs)
            peaks_dict[n] = _peakslocs
        has_peak_in_utr = False
        for _peak in txpts[gene_name].peak_locs:
            for pos in range(_peak[0], _peak[1] + 1):
                pos_in_right = pos - txpts[gene_name].stop_pos_in_seq
                if 0 < pos_in_right < len(_utrsarr):
                    has_peak_in_utr = True
                    _peaks_in_utr[pos_in_right] = np.log10(_peak[2])
        if has_peak_in_utr:
            utrs_dict[n] = _utrsarr
            peaks_in_utrs_dict[n] = _peaks_in_utr
    return (txpts_dict, peaks_dict, utrs_dict, peaks_in_utrs_dict)


def has_peak_in_range(txpts, gene_name, max_value):
    for _peak in txpts[gene_name].peak_locs:
        for pos in range(_peak[0], _peak[1] + 1):
            if(0 < pos < max_value):
                return True
    return False


def add_peak_signal(_txpts, gene_name, _peakslocs):
    for _peak in _txpts[gene_name].peak_locs:
        for pos in range(_peak[0], _peak[1] + 1):
            if(0 < pos < len(_peakslocs)):
                _peakslocs[pos] = np.log10(_peak[2])
