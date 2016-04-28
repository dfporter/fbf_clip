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


def sort_and_make_fractions(peaks):
    sorted_peaks = sorted(peaks.items(), key=operator.itemgetter(0))
    sum_peaks = float(sum([x[1] for x in sorted_peaks]))
    sorted_peaks = [(x[0], float(x[1])/sum_peaks) for x in sorted_peaks]
    return sorted_peaks


def write_raw_table(list_of_lists):
    li = 'Position\tMotif\n'
    for row in list_of_lists:
        li += '{a}\t{b}\n'.format(
            a=row[0], b=row[1]
        )
    return li


def plot_features(aves, output_filename='../clip/figs/features_in_normalized_gene.pdf'):
    (ave_peaks, ave_fbes, ave_negs, ave_polyA, ave_highest_peak, ave_secondary_peaks) = aves
    plt.clf()
    sorted_peaks = sort_and_make_fractions(ave_peaks)
    sorted_highest_peak = sort_and_make_fractions(ave_highest_peak)
    sorted_secondary_peaks = sort_and_make_fractions(ave_secondary_peaks)
    sorted_fbes = sort_and_make_fractions(ave_fbes)
    sorted_negs = []
    for neg in ave_negs:
        #print 'init: %s' % str(neg)
        ret = sort_and_make_fractions(neg)
        #print 'ret: %s' % str(ret)
        sorted_negs.append(ret)
#    sys.exit()
    sorted_polyA = sort_and_make_fractions(ave_polyA)
    plt.xlabel('Position (nt)')
    plt.ylabel('Frequency')
    plt.plot([x[0] for x in sorted_fbes],
                [x[1] for x in sorted_fbes], label='FBE')
#    for index, neg in enumerate(sorted_negs):
#        plt.plot([x[0] for x in neg],
#                 [x[1] for x in neg], linestyle='--', label='neg-' + str(index))
#    plt.plot([x[0] for x in sorted_polyA],
#                [x[1] for x in sorted_polyA], label='poly(A)')
    plt.plot([x[0] for x in sorted_peaks],
                [x[1] for x in sorted_peaks], label='Peaks')
    with open('peaks_in_normalized_gene.txt', 'w') as f:
        f.write(write_raw_table(sorted_peaks))
    with open('motifs_in_normalized_gene.txt', 'w') as f:
        f.write(write_raw_table(sorted_fbes))

    pos = np.arange(0, 1400 + 200,200)
    _ylabels = [r'''Start of UTR''', r'''Start''',
                '', '', '', '',
                r'''Stop''', r'''End of UTR''']
    plt.xticks(pos, _ylabels)
    plt.legend(loc='upper left')
    print "\tSaving figure to %s" % output_filename
    plt.savefig(output_filename, format='pdf')
    plt.clf()
    # Zoom in on UTR.
    utr_output_filename = output_filename.partition('.pdf')[0] + '_zoom.pdf'
    plt.xlabel('Position (nt)')
    plt.ylabel('Frequency')
    plt.plot([x[0] for x in sorted_fbes],
                [x[1] for x in sorted_fbes], label='FBE')
    plt.plot([x[0] for x in sorted_peaks],
                [x[1] for x in sorted_peaks], label='Peaks')
    plt.xlim([1200,1400])
    pos = np.arange(1200, 1400 + 100, 100)
    _ylabels = [r'''Start of UTR''', '',  r'''End of UTR''']
    plt.xticks(pos, _ylabels)
    plt.legend(loc='upper left')
    print "\tSaving figure to %s" % output_filename
    plt.savefig(utr_output_filename, format='pdf')
    with open('peaks_in_normalized_utr_zoom.txt', 'w') as f:
        f.write(write_raw_table(sorted_peaks[1200:]))
    plt.clf()
    # Highest peak vs secondary peaks.
    plt.xlabel('Position (nt)')
    plt.ylabel('Frequency')
    plt.plot([x[0] for x in sorted_highest_peak],
                [x[1] for x in sorted_highest_peak], label='Highest peak')
#    plt.plot([x[0] for x in sorted_polyA],
#                [x[1] for x in sorted_polyA], label='poly(A)')
    plt.plot([x[0] for x in sorted_secondary_peaks],
                [x[1] for x in sorted_secondary_peaks], label='Secondary peaks')
    plt.legend(loc='upper left')
    pos = np.arange(0, 1400 + 200,200)
    _ylabels = [r'''Start of UTR''', r'''Start''',
                '', '', '', '',
                r'''Stop''', r'''End of UTR''']
    plt.xticks(pos, _ylabels)
    output_filename_secondary = output_filename.partition('.pdf')[0] + '_secondary.pdf'
    plt.savefig(output_filename_secondary, format='pdf')
    plt.clf()


def plot_utr(ave_peak, ave_fbe, output_filename='../clip/figs/features_in_utr.pdf'):
    sorted_peaks = sort_and_make_frequency(ave_peak)
    sorted_fbe = sort_and_make_frequency(ave_fbe)
    plt.clf()
    # Peaks.
    plt.xlabel('Position (nt)')
    plt.ylabel('Frequency')
    plt.plot([x[0] for x in sorted_peaks],
                [x[1] for x in sorted_peaks],'b-', label='Peaks')
    plt.plot([x[0] for x in sorted_fbe],
                [x[1] for x in sorted_fbe], 'r-', label='FBE')
    with open('peaks_in_utr_raw_data.txt', 'w') as f:
        f.write(write_raw_table(sorted_peaks))
    with open('motifs_in_utr_raw_data.txt', 'w') as f:
        f.write(write_raw_table(sorted_fbe))
    plt.legend(loc='upper left')
    plt.xticks(range(0,350,50), ['Stop', '-150', '-100', '-50', 'AAUAAA', '50', 'End of UTR'])
    plt.savefig(output_filename, format='pdf')
    plt.clf()


def sort_and_make_frequency(aves):
    sorted_peaks = sorted(aves.items(), key=operator.itemgetter(0))
    sum_peaks = float(sum([x[1] for x in sorted_peaks]))
    sorted_peaks = [(x[0], float(x[1])/sum_peaks) for x in sorted_peaks]
    return sorted_peaks


def normalize_distances(txpts):
    for t in txpts:
        print 'norm:'
        print t
        txpts[t].normalize_features()
        txpts[t].get_peak_pos_relative_to_polyA()
        txpts[t].get_motif_pos_relative_to_polyA()


def get_ave_pos(txpts):
    if len(txpts) == 0:
        print "get_ave_pos(): Error, empty txpts."
        return
    ave_motif = {}
    ave_polyA = {}
    ave_peaks = {}
    ave_highest_peak = {}
    ave_secondary_peak = {}
    ave_negs = [{}, {}, {}]
    # FBEs, polyA, and other motifs.
    print "*(*(*" * 17
    for t in txpts:
        print t
        print txpts[t].motif_locs
        print txpts[t].norm_motif
        for _span in txpts[t].norm_motif:
            print 'asdf'
            print _span
            for pos in range(_span[0], _span[1]):
                ave_motif.setdefault(pos, 0.)
                ave_motif[pos] += 1.
        for index, neg in enumerate([txpts[t].norm_neg1, txpts[t].norm_neg2, txpts[t].norm_neg3]):
            for _span in neg:
                for pos in range(_span[0], _span[1]):
                    ave_negs[index].setdefault(pos, 0.)
                    ave_negs[index][pos] += 1
        for _span in txpts[t].norm_polyA:
            print 'a'
            print _span
            for pos in range(_span[0], _span[1]):
                ave_polyA.setdefault(pos, 0.)
                ave_polyA[pos] += 1.
    for pos in ave_motif:
        ave_motif[pos] = ave_motif[pos]/float(len(txpts))
    for neg in ave_negs:
        for pos in neg:
            neg[pos] = neg[pos]/float(len(txpts))
    for pos in ave_polyA:
        ave_polyA[pos] = ave_polyA[pos]/float(len(txpts))
    # Peaks.
    for t in txpts:
        for _span in txpts[t].norm_peaks:
            for pos in range(_span[0], _span[1]):
                ave_peaks.setdefault(pos, 0.)
                ave_peaks[pos] += float(_span[2])  # Height.
        if len(txpts[t].norm_peaks) > 1:
            sorted_peaks = sorted(txpts[t].norm_peaks, key=lambda x: float(x[2]))
            highest_peak = sorted_peaks[-1]
            secondary_peaks = sorted_peaks[0:-1]
            for pos in range(highest_peak[0], highest_peak[1]):
                for pos in range(highest_peak[0], highest_peak[1]):
                    ave_highest_peak.setdefault(pos, 0.)
                    ave_highest_peak[pos] += highest_peak[2]
            for _span in list(secondary_peaks):
                for pos in range(_span[0], _span[1]):
                    ave_secondary_peak.setdefault(pos, 0.)
                    ave_secondary_peak[pos] += _span[2]
    for pos in ave_peaks:
        ave_peaks[pos] = ave_peaks[pos]/float(len(txpts))
    for pos in ave_highest_peak:
        ave_highest_peak[pos] = ave_highest_peak[pos]/float(len(txpts))
    for pos in ave_secondary_peak:
        ave_secondary_peak[pos] = ave_secondary_peak[pos]/float(len(txpts))
    ave_peak_v_polyA = {}
    for t in txpts:
        if not hasattr(txpts[t], 'norm_peak_v_polyA'): continue
        if not txpts[t].norm_peak_v_polyA: continue
        for pos in range(int(txpts[t].norm_peak_v_polyA[0]),
                         int(txpts[t].norm_peak_v_polyA[1])):
            ave_peak_v_polyA.setdefault(pos, 0.)
            ave_peak_v_polyA[pos] += float(txpts[t].norm_peak_v_polyA[2])
    ave_motif_v_polyA = {}
    for t in txpts:
        if not hasattr(txpts[t], 'norm_motif_v_polyA'):
            continue
        if not txpts[t].norm_motif_v_polyA:
            continue
        for _norm_p in txpts[t].norm_motif_v_polyA:
            if not _norm_p: continue
            for pos in range(int(_norm_p[0]),
                             int(_norm_p[1])):
                ave_motif_v_polyA.setdefault(pos, 0.)
                ave_motif_v_polyA[pos] += 1.0
    return [ave_peaks, ave_motif, ave_negs, ave_polyA, ave_peak_v_polyA,
            ave_motif_v_polyA,
            ave_highest_peak, ave_secondary_peak]

def peak_vs_motif(txpts):
    ave_dist = {}
    for t in txpts:
        if (len(txpts[t].norm_peaks) > 0) and (len(txpts[t].norm_motif) > 0):
            for dist in get_closest_dist(
                    txpts[t].norm_peaks, txpts[t].norm_motif):
                ave_dist.setdefault(dist, 0.)
                ave_dist[dist] += 1.


def get_closest_dist(tups1, tups2):
    closest = []
    for tup in tups1:
        closest.append(100000000.)
        for tup2 in tups2:
            nearest = min(abs(tup[0] - tup2[0]),
                abs(tup[0] - tup2[1]),
                abs(tup[1] - tup2[0]),
                abs(tup[1] - tup2[1]))
            if nearest < closest[-1]:
                closest[-1] = nearest
    return closest

