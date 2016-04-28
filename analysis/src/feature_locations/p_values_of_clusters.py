"""

0. For a give target, get longest txpt.
1. Get all exon genomic ranges for txpt.
2. Determine clusters (start, stop, max overlap)
3. For the total number of reads, scramble their positions within the exon ranges.
4. Get a histogram from the scrambles.

"""
import HTSeq
import numpy as np
import bisect
import pandas
# from collections import Counter
import collections
import copy
import os


def generate_ivs(num_reads=100, start=0, end=40, chrom='I', strand='+', ranges=None,
                 fast=False):
    """
    If fast=True, don't return HTSeq.GenomicIntervals.
    """
    ivs = []
    if ranges is not None:
        total_len_of_ranges = 0
        list_of_range_starts = []
        ranges = sorted(ranges, key=lambda x: x[0])
        left_edge = ranges[0][0]
        for _range in ranges:
            list_of_range_starts.append(_range[0])
            total_len_of_ranges += _range[1] - _range[0]
        if fast:
            reads = []
            for n in np.random.random_integers(0, total_len_of_ranges, size=num_reads):
                reads.append([n, n + np.random_integers(5, high=20)])
            return reads
        for i in range(num_reads):
            rand_start = np.random.random_integers(
                left_edge, high=left_edge+total_len_of_ranges)
            start_exon_index = bisect.bisect(list_of_range_starts, rand_start) - 1
            if start_exon_index < 0:
                continue
            start_rel_pos = rand_start + list_of_range_starts[start_exon_index]
            end_rel_pos = start_rel_pos + \
                np.random.random_integers(5, high=20)
            ivs.append(HTSeq.GenomicInterval(chrom, start_rel_pos, end_rel_pos, strand))
        return ivs, total_len_of_ranges
    # for i in range(n):
    #     rand_start = np.random.random_integers(start, high=end)
    #     rand_end = rand_start + \
    #         np.random.random_integers(5, high=20)
    #     ivs.append(HTSeq.GenomicInterval(chrom, rand_start, rand_end, strand))
    return ivs


def get_ranges_set_by_coverage_depth(points):
    ranges = []  # Holds (start, stop, set of indexes) for output.
    current_set = set()
    last_start = None
    for position, cat, index in points:
        if cat == 'start':
            if last_start is not None:
                ranges.append((last_start, position, copy.copy(current_set)))
            current_set.add(index)
            last_start = position
        elif cat == 'end':
            ranges.append((last_start, position, copy.copy(current_set)))
            current_set.remove(index)
            last_start = position
    return ranges


def fast_find_clusters_in_gene(positions):
    clusters = []
    current_cluster = None
    max_overlap_in_cluster = 0
    positions = sorted(positions, key=lambda x: x[0])
    for a_range in positions:
        if current_cluster is None:
            current_cluster = [a_range[0], a_range[1]]
            max_overlap_in_cluster = 1
        else:
            # Do we overlap the current cluster?
            if a_range[0] < current_cluster[1]:
                max_overlap_in_cluster += 1
                current_cluster[1] = a_range[1]
            else:
                clusters.append(current_cluster + [max_overlap_in_cluster])
                current_cluster = [a_range[0], a_range[1]]
                max_overlap_in_cluster = 1
    if current_cluster is not None:
        clusters.append(current_cluster + [max_overlap_in_cluster])
    if len(clusters) == 0:
        clusters = [[0, 0, 1]]
    return clusters


def find_clusters_in_gene(rigs):
    """
    Find overlapping reads in a given gene from a list of ivs.
Input:
    rigs: A list of iv objects, one for each read. Includes all reads in a given
            gene.
Internal:
    points: A list of (start, stop, read index) for each start and stop,
            made separately for each read.
    ranges: A list of (start, stop, set of read indexes in range).
Output:
    clusters: A list of [start, stop, max read overlap] for each region of
            continuous coverage.
    """
    # rigs is a list of HTSeq iv objects, in genomic coordinates.
    points = []  # Holds (start, stop, symbol).
    for index, iv in enumerate(rigs):
        points.append((iv.start, 'start', index))
        points.append((iv.end, 'end', index))
    points.sort()
    ranges = get_ranges_set_by_coverage_depth(points)
    # ranges = [(start, stop, set of indexes)..].
    clusters = []
    current_cluster = None
    max_overlap_in_cluster = 0
    for a_range in ranges:
        if current_cluster is None:
            if len(a_range[2]) > 0:
                current_cluster = [a_range[0], a_range[1], len(a_range[2])]
                max_overlap_in_cluster = len(a_range[2])
        else:
            if len(a_range[2]) > 0:
                current_cluster[1] = a_range[1]
                max_overlap_in_cluster = max([max_overlap_in_cluster, len(a_range[2])])
                current_cluster[2] = max_overlap_in_cluster
            else:
                max_overlap_in_cluster = max([max_overlap_in_cluster, len(a_range[2])])
                current_cluster[1] = a_range[0] - 1
                current_cluster[2] = max_overlap_in_cluster
                clusters.append(current_cluster)
                current_cluster = None
    if current_cluster is not None:
        current_cluster[2] = max_overlap_in_cluster
        clusters.append(current_cluster)
    # print "*"
    # print clusters
    # if last_start is not None:
    #     ranges.append((last_start, position, copy.copy(current_set)))
    return clusters
#    visualize(rigs, ranges, clusters)


def visualize(rigs, ranges, clusters):
    li = ''
    for y in range(0, 6):
        for x in range(0, 10):
            if x == 9:
                li += str(y)
            else:
                li += ' '
    print li
    print "".join([str(x) for x in range(10)]) * 6
    for clust in clusters:
        print " " * clust[0] + "c" * (clust[1] - clust[0])
    for iv in rigs:
        print str(
            " " * iv.start + "#" * (iv.end - iv.start)
        )
    for tup in ranges:
        print str(
            " " * tup[0] + "=" * (tup[1] - tup[0])
        ) + str([str(x) for x in tup[2]])


def get_exonic_ranges(fname):
    import sys

    df = pandas.read_csv(fname, sep='\t')
    gtfd = df.to_dict('records')
    by_gene = collections.defaultdict(list)
    for row in gtfd:
        if row['2'] != 'exon':
            continue
        by_gene[row['gene_id']].append(row)
    return by_gene


def get_exon_ranges_and_scramble_and_return_probabilities(
        by_gene=None, gene_name=None, num_reads=2, num_permutations=1000,
        txpt_id=None, exons=None):
    if txpt_id is None or exons is None:
        txpt_id, exons = find_longest_txpt(by_gene[gene_name])  # Returns txpt_id, [(start, stop)...]
    max_coverages = []
    total_len_of_ranges = 0
    list_of_range_starts = []
    ranges = sorted(exons, key=lambda x: x[0])
    for _range in ranges:
        list_of_range_starts.append(_range[0])
        total_len_of_ranges += _range[1] - _range[0]
    if num_reads>1:
        for j in range(num_permutations):
            reads = []
            rand_starts = np.random.random_integers(0, total_len_of_ranges, size=num_reads)
            for n in rand_starts:
                reads.append([n, n + np.random.random_integers(5, high=20)])
            clusters = fast_find_clusters_in_gene(reads)
            max_coverages.append(max([x[2] for x in clusters]))
        return get_histogram(max_coverages)
    for j in range(num_permutations):
        ivs, total_len = generate_ivs(
            num_reads=num_reads, ranges=exons)
        clusters=find_clusters_in_gene(ivs)
        if len(clusters) > 0:
            max_coverages.append(max([x[2] for x in clusters]))
        else:
            max_coverages.append(1)
    return get_histogram(max_coverages)


def get_histogram(observed_depths, verbose=False):
    # _counts = Counter(observed_depths)  # Need python2.7+
    _counts = collections.defaultdict(int)
    for obs in observed_depths:
        _counts[obs] += 1
    keys = sorted(_counts.keys())
    total_tests = float(len(observed_depths))
    sum_so_far = 0.
    probs = collections.defaultdict(float)
    for key in keys:
        probs[key] = float(
            1 - float(sum_so_far/max([total_tests, 1.]))
        )
        if verbose:
            print """\
{k}: {i} (cumulative: {s}/{t}) -> probability of this value or higher {p}""".format(
            k=key, i=_counts[key], s=sum_so_far, t=total_tests,
            p=probs[key],
        )
        sum_so_far += float(_counts[key])
    return probs


def find_longest_txpt(exons, gene_id='Unknown'):
    longest_txpt = None
    known_txpts = collections.defaultdict(list)
    for row in exons:
        known_txpts[row['transcript_id']].append([row['3'], row['4']])
    txpt_lens = collections.defaultdict(int)
    for txpt in known_txpts:
        txpt_lens[txpt] = 0
        for exon in known_txpts[txpt]:
            txpt_lens[txpt] += int(exon[1]) - int(exon[0])
    longest_txpt = sorted(txpt_lens.keys(), key=lambda x: txpt_lens[x])[-1]
    li = "Txpts for {a}\n".format(a=gene_id)
    for key in known_txpts:
        li += """{i}: {le}\n""".format(i=key, le=txpt_lens[key])
    li += 'Picked {k}\n'.format(k=longest_txpt)
    # print li
    return longest_txpt, known_txpts[longest_txpt]


def add_p_values_from_histogram(probabilities, clusters):
    """
    Find overlapping reads in a given gene from a list of ivs.
Input:
    clusters: A list of [start, stop, max read overlap] for each region of
            continuous coverage in a given gene.
    probabilites: (This is for a given transcript). A dict mapping a given
            count of maximum read coverage to the probability of obtaining
            a cluster at least that high. This is obtained from permutations,
            so not every count is obtained.
Output:
    clusters_w_pvalue: A list of [start, stop, max read overlap, p value].
            It's the same as the input clusters list, but with a p value.
    """
    clusters_w_pvalue = []
    counts_with_known_p_values = sorted(probabilities.keys(), key=lambda x: x)
    for clust in clusters:
        clusters_w_pvalue.append(clust)
        index = bisect.bisect(counts_with_known_p_values, clust[2])
        if index == 0:
            clusters_w_pvalue[-1].append(1.)
        elif index >= len(counts_with_known_p_values):
            clusters_w_pvalue[-1].append(0.)
        else:
            closest_measured_value = counts_with_known_p_values[index-1]
            clusters_w_pvalue[-1].append(probabilities[closest_measured_value])
    return clusters_w_pvalue