import re
import sys
import pandas
import numpy as np
import datetime
import logging

#fname = sys.argv[1]

# Want to compare single peaks vs peak pairs for their locations,
# specifically whether the tendency to be in the 3'UTR is different between the two.
# Second - what proportion of peaks (primary, secondary, both in the case of a pair).

# So we want a list of the locations for each category.
# One could build the 2x2 contingency table:
#       is primary
#       -       +
#    -
#  in 3'UTR
#    +
# And the same table for peak pairs vs single peaks.



def count_loc(subset_of_peaks):
    in_utr = 0
    not_in_utr = 0
    for gene in subset_of_peaks:
        for peak in subset_of_peaks[gene]:
            if peak['location'] == "3'UTR":
                in_utr += 1
            else:
                not_in_utr += 1
    return(in_utr, not_in_utr)

def convert_df_to_dict(df):
    peaks_d = {}
    for index, row in df.iterrows():
        peaks_d.setdefault(row['gene_name'], [])
        peaks_d[row['gene_name']].append(row.to_dict())
    return peaks_d


def get_gene_size(gtf):
    gene_lens = {}
    for n, gene_name in enumerate(gtf):
        # Find the longest txpt.
        txpt_ids = set([x['transcript_id'] for x in gtf[gene_name]]) \
                   - set([np.nan, '.'])
        txpts_exons = {}
        for t in txpt_ids:
            txpts_exons.setdefault(
                t, [x for x in gtf[gene_name] if (
                    (x['transcript_id']==t) and (x['2']=='exon'))]
            )
        if len(txpts_exons) == 0: continue
        longest_txpt_id = ('na', -1)
        for t in txpts_exons:
            t_right = max([int(x['4']) for x in txpts_exons[t]])
            t_left = min([int(x['3']) for x in txpts_exons[t]])
            t_len = t_right - t_left
            if t_left > t_right:
                print "Coordinate error in %s" % str(gtf[gene_name])
                return False
            if t_len > longest_txpt_id[1]:
                longest_txpt_id = (t, t_len)
        # Measure gene feature lengths in the longest transcript.
        exons = txpts_exons[longest_txpt_id[0]]
        cds = [x for x in gtf[gene_name] if (
            (x['transcript_id'] == longest_txpt_id[0]) and (x['2']=='CDS'))]
        if len(cds) == 0: continue
        min_cds = min([int(x['3']) for x in cds])
        max_cds = max([int(x['4']) for x in cds])
        min_exon = min([int(x['3']) for x in exons])
        max_exon = max([int(x['4']) for x in exons])
        left_utr = min_cds - min_exon
        right_utr = max_cds - max_exon
        gene_lens[gene_name] = {}
        if gtf[gene_name][0]['6'] == '+':
            gene_lens[gene_name]["5'UTR"] = left_utr
            gene_lens[gene_name]["3'UTR"] = right_utr
        if gtf[gene_name][0]['6'] == '-':
            gene_lens[gene_name]["3'UTR"] = left_utr
            gene_lens[gene_name]["5'UTR"] = right_utr
        gene_lens[gene_name]['txpt_len'] = 0
        for a_exon in exons:
            gene_lens[gene_name]['txpt_len'] += int(a_exon['4']) - int(a_exon['3'])
        gene_lens[gene_name]['CDS'] = 0
        for a_exon in cds:
            gene_lens[gene_name]['CDS'] += int(a_exon['4']) - int(a_exon['3'])
    print "Lengths for %i genes loaded." % len(gene_lens)
#    print "Examples: %s" % str([
#        gene_lens[gene_lens.keys()[x]] for x in range(0, min([len(gene_lens), 10]))
#        ])
    return gene_lens



def write_line(tup):
    return "In 3'UTR: {i} Outside 3'UTR: {o} (sum={s})\n\t(Percent in 3'UTR: {p})".format(
        i=tup[0], o=tup[1], s=sum(tup), p="%.2f" % (100. * float(tup[0])/float(sum(tup)))
    )

def read_peaks_file(fname):
    peaks = pandas.read_csv(fname, sep='\t')
    if 'location' not in peaks.columns:
        print "Did not find location column in %s. Adding it..." % fname
        gtf = pandas.read_csv(fname, sep='\t')
        import add_info_columns
        add_info_columns.locate_in_gene(gtf, peaks)
        peaks.to_csv(fname, sep='\t', index=False)
    by_gene = {}
    known_indexes = set()
    peaks_read = 0
    for index, row in peaks.iterrows():
        if index in known_indexes:
            print "error: duplicate indexes"
            return False
        by_gene.setdefault(row['gene_name'], [])
        by_gene[row['gene_name']].append(row.to_dict())
        peaks_read += 1
    print "Total: {i} peaks...".format(i=peaks_read)
    return by_gene


def run(peaks_filename, gtf_d):
    #peaks_filename = '/groups/Kimble/Aman Prasad/clip/methods/filtered_gauss_same_n2_nil_five_reps/combined_fbf.txt'
    #peaks_filename = '/groups/Kimble/Aman Prasad/clip/methods/filtered_gauss_nil_five_reps/combined_fbf.txt'

    gene_lens = get_gene_size(gtf_d)
    #(txpt_size, utr_size) = get_gene_size(gtf_d)
    peaks_by_gene = read_peaks_file(peaks_filename)
    where_are_peak_pairs(peaks_by_gene)
    size_distrubtion_of_txpt_len_for_peak_type(peaks_by_gene, gene_lens)
    tups = separate_genes_by_peak_number(peaks_by_gene)
    labels = ('single_peak', 'primary_peak', 'peak_pairs', 'secondary_peak')
    split_peaks = {}
    import scipy.stats
    for n, k in enumerate(tups):
        split_peaks[labels[n]] = k
    arrs = {}
    for k in split_peaks:
        print k
        arrs[k] = size_distrubtion_of_txpt_len_for_peak_type(split_peaks[k], gene_lens)
    print 't-test for txpt length differences between genes with single peaks or multiple peaks'
    print scipy.stats.ttest_ind(arrs['single_peak']['txpt_len'],
                                arrs['primary_peak']['txpt_len'])
    print 't-test for cds length differences between genes with single peaks or multiple peaks'
    print scipy.stats.ttest_ind(arrs['single_peak']['CDS'],
                                arrs['primary_peak']['CDS'])

    #num_pks = count_all_locs(single=single_peak, primary_peak=primary_peak,
    #               peak_pairs=peak_pairs, seconday_peak=secondary_peak)


def size_distrubtion_of_txpt_len_for_peak_type(peaks_by_gene, gene_len):
    #for gene in by_gene:
    #    by_gene[gene].update({
    #        'txpt_len': gene_len[gene]['txpt_len'],
    #        '5utr_size': gene_len[gene]["5'UTR"],
    #        '3utr_size': gene_len[gene]["3'UTR"],
    #        'cds_size': gene_len[gene]["CDS"],
    #        })
    darrs = {"3'UTR": {}, "5'UTR": {}, 'CDS': {}, 'txpt_len': {}}
    for gene in peaks_by_gene:
        if gene not in gene_len:
            #print 'missing %s' % gene
            continue
        for peak in peaks_by_gene[gene]:
            if peak['location'] == """3'UTR""":
                if "3'UTR" in gene_len[gene]:
                    darrs["3'UTR"][gene] = gene_len[gene]["3'UTR"]
                else: "no UTR"
            if peak['location'] == "CDS":
                if "CDS" in gene_len[gene]:
                    darrs["CDS"][gene] = gene_len[gene]["CDS"]
                else: 'no UTR'
            if "txpt_len" in gene_len[gene]:
                darrs["txpt_len"][gene] = gene_len[gene]["txpt_len"]
            else: 'no txpt_len'
    arrs = {}
    arrs["3'UTR"] = sorted(darrs["3'UTR"].values())
    arrs["CDS"] = sorted(darrs["CDS"].values())
    arrs["txpt_len"] = sorted(darrs["txpt_len"].values())
    print "\tNumbers of features (# genes) counted: %s" % str(
        (arrs.keys(), [len(arrs[k]) for k in arrs]))
    for k in arrs:
        print "\t{k}: mean {mn} median {md}".format(
            k=k, mn=np.mean(arrs[k]), md=np.median(arrs[k]))
    return arrs


def separate_genes_by_peak_number(by_gene):
    single_peak = {}
    primary_peak = {}
    peak_pairs = {}
    secondary_peak = {}
    for gene_name in by_gene:
        if len(by_gene[gene_name]) == 1:
            single_peak[gene_name] = [by_gene[gene_name][0]]
        elif len(by_gene[gene_name]) > 1:
            sortedpeaks = sorted(
                by_gene[gene_name],
                key=lambda x: int(x['height']))
            primary_peak[gene_name] = [sortedpeaks[-1]]
            for p in sortedpeaks[:-1]:
                secondary_peak.setdefault(gene_name, [])
                secondary_peak[gene_name].append(p)
            peak_pairs[gene_name] = sortedpeaks
    return single_peak, primary_peak, peak_pairs, secondary_peak


def where_are_peak_pairs(by_gene):
    assert type(by_gene) is type({})
    if len(by_gene) > 0: assert type(by_gene.values()[0]) is type([])
    single_peak = {}
    peak_pairs = {}
    for gene_name in by_gene:
        if len(by_gene[gene_name]) == 1:
            single_peak[gene_name] = [by_gene[gene_name][0]['location']]
        elif len(by_gene[gene_name]) > 1:
            sortedpeaks = sorted(
                by_gene[gene_name],
                key=lambda x: int(x['height']))
            peak_pairs[gene_name] = []
            for p in sortedpeaks:
                peak_pairs[gene_name].append(p['location'])
    all_in_3utr = set()
    not_all_in_3utr = set()
    for gene in peak_pairs:
        if set(peak_pairs[gene]) == set(["3'UTR"]):
            all_in_3utr |= set([gene])
        else:
            not_all_in_3utr |= set([gene])
    fraction = float(len(list(all_in_3utr)))/float(len(peak_pairs))
    print """
    Peak pairs {a}
    Peak pairs all in 3'UTR: {u}
    Peak pairs all in 3'UTR: {f}""".format(
        a=len(peak_pairs),
        u=len(list(all_in_3utr)),
        f=fraction)


def count_all_locs(**kwargs):
    print "Counting peak locations."
    num_pks = {}
    for key in kwargs:
        num_pks[key] = count_loc(kwargs[key])
        print key + ': ' + write_line(num_pks[key])
    return num_pks
