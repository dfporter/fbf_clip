import pandas
import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys
import glob
import numpy as np
import re
import scipy
import scipy.stats


import compare_with_ripchip

def get_height(row, ga):
    left = row['left']
    right = row['right']
    chrm = row['chrm']
    strand = row['strand']
    iv = HTSeq.GenomicInterval(
        chrm, left, right, strand
    )
    coverage = np.fromiter(ga[iv], dtype='i')
    height = np.max(coverage)
    return height


def for_every_row_get_heights_over_motifs(peaks, ga, motifs):
    peaks_d = peaks.to_dict('records')
    signals = {}
    chr_lens = {}
    with open('lib/chr_sizes.txt', 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            chr_lens[s[0]] = int(s[1])
    for row in peaks_d:
        get_height_over_motifs(row, motifs, chr_lens, ga, signals=signals)
    pos_arrs = {}
    for key in signals:
        positives = filter(lambda x: len(x) > 0, signals[key])
        pos_arr = []
        for row in positives:
            for num in row:
                if num > 0:
                    pos_arr.append(num)
        print key
        #print pos_arr
        print np.median(pos_arr)
        print np.mean(pos_arr)
        pos_arrs[key] = pos_arr
    print scipy.stats.ttest_ind(pos_arrs['tgt\w\w\wat'], pos_arrs['tgt[ag]aat'])


def get_height_over_motifs(row, motifs, chr_lens, ga, signals={}):
    left = row['left']
    right = row['right']
    chrm = row['chrm']
    strand = row['strand']
    #iv = HTSeq.GenomicInterval(
    #    chrm, left, right, strand
    #)
    chr_len = chr_lens[chrm]
    for motif in motifs: signals.setdefault(motif, [])
    for motif in motifs:
        signals[motif].append(
            get_height_over_a_motif(motif, ga, row, chr_len)
        )
    return signals

def compare_heights_over_motifs(peaks, ga, motifs):
    peaks_d = peaks.to_dict('records')
    signals = {}
    chr_lens = {}
    with open('lib/chr_sizes.txt', 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            chr_lens[s[0]] = int(s[1])
    for row in peaks_d:
        get_height_over_motifs(row, motifs, chr_lens, ga, signals=signals)
    pos_arrs = {}
    for key in signals:
        positives = filter(lambda x: len(x) > 0, signals[key])
        pos_arr = []
        for row in positives:
            for num in row:
                if num > 0:
                    pos_arr.append(num)
        print key
        #print pos_arr
        print np.median(pos_arr)
        print np.mean(pos_arr)
        pos_arrs[key] = pos_arr
    print scipy.stats.ttest_ind(pos_arrs['tgt\w\w\wat'], pos_arrs['tgt[ag]aat'])

def get_height_over_a_motif(motif, ga, row, chr_len):
    motif_pos = []
    motif_pos_genomic = []
    for m in re.finditer(motif, row['seq']):
            motif_pos.append((m.start(), m.end()))
    for pos in motif_pos:
        if row['strand'] == '+':
            motif_pos_genomic.append(
                HTSeq.GenomicInterval(row['chrm'], row['left'] + pos[0], row['left'] + pos[1],
                                      row['strand'])
            )
        if row['strand'] == '-':
            motif_pos_genomic.append(
                HTSeq.GenomicInterval(row['chrm'], row['right'] - pos[1],
                                      row['right'] - pos[0],
                                      row['strand'])
            )
    max_signals = []
    for iv in motif_pos_genomic:
        max_signals.append(
            (np.max(np.fromiter(ga[iv], dtype='i'))
             )
        )
    return max_signals

def sort_cols(df):
    header = [
        'chrm', 'left', 'right', 'name', 'gene_name', 'strand',
        'no_antibody_ratio', 'n2_ratio', 'ratio', 'gene_id', 'biotype',
        'location', 'RIP_target', 'RIP_rank', 'height', 'padj', 'pvalue',
        'log2FoldChange', 'seq', 'fog_CGGA', 'fog_GGCA', 'fog_GGTT', 'fog_TGGC',
        'control_AATA', 'control_CCGG', 'control_TTAA', 'control_TTGT',
        'fbf1n2_acaa', 'fbf1n2_ccgg', 'fbf1n2_tgcc', 'fbf2n2_acaa',
        'fbf2n2_ccgg', 'fbf2n2_tgcc']
    unrecognized = list(set(df.columns) - set(header))
    sort_by = []
    for col in header:
        if col in df.columns:
            sort_by.append(col)
    sort_by.extend(unrecognized)
    df = df[sort_by]
    return df


def score_positives(peaks):
    known_pos = set(['gld-1', 'htp-1', 'htp-2', 'mpk-1', 'him-3',
                         'fbf-1', 'lip-1', 'syp-2', 'fbf-2', 'fog-1',
                         'fem-3', 'syp-3', 'gld-3', 'fog-3', 'egl-4'])
    if 'Rank' not in peaks.columns:
        peaks.sort(columns='fbf2_reads', ascending=False, inplace=True)
        peaks['Rank'] = range(len(peaks.index))
    obs_genes = set(peaks['gene_name'].tolist())
    obs_pos = known_pos & obs_genes
    missing_pos = known_pos - obs_genes
    obs_pos_n = len(list(obs_pos))
    missing_pos_n = len(list(missing_pos))
    best_ranks = {}
    for gene in obs_pos:
        rows = peaks[peaks['gene_name']==gene]
        rows = rows.to_dict('records')
        best_rank = 1e4
        for a_row in rows:
            best_rank = min([a_row['Rank'], best_rank])
        best_ranks[gene] = best_rank
    num_below_500 = 0
    below_500 = filter(lambda x: x<=500, best_ranks.values())
    print below_500
    below_700 = filter(lambda x: x<=700, best_ranks.values())
    print below_700
    print scipy.stats.fisher_exact([[obs_pos_n, missing_pos_n],
                                    [len(list(obs_genes)), 2e4 - missing_pos_n - obs_pos_n - len((list(obs_genes)))]
                                    ])
    return {'observed positives': obs_pos, 'number of observed positives': obs_pos_n,
            'missing positives': missing_pos, 'number of missing positives': missing_pos_n,
            'expected': len(list(known_pos)),
            'ranks:': best_ranks, 'Positives below rank 700:': len(below_700)}


def score_fbe(peaks):
    peaks = peaks[peaks['biotype']=='protein_coding']
    with_fbe = peaks[peaks['has_fbe']==1]
    without_fbe = peaks[peaks['has_fbe']==0]
    frac_with = float(len(with_fbe.index))/float(
        len(with_fbe.index) + len(without_fbe.index))
    frac_without = float(len(without_fbe.index))/float(
        len(with_fbe.index) + len(without_fbe.index))
    print "Peaks with FBE: {w} ({wp}%). Peaks without FBE: {wo} ({a}%).".format(
        w=len(with_fbe.index), wo=len(without_fbe.index),
        wp=100*frac_with, a=100*frac_without
    )


def split_fasta(df, filename='test.fa'):
    li = ''
    for index, row in df.iterrows():
        li += '>{i}_{n}\n'.format(row['Rank'], row['gene_name'])
        li += ''

def write_fasta_by_location(peaks):
    seqs = {}
    stats = {}
    for index, row in peaks.iterrows():
        seqs.setdefault(row['location'], '')
        seqs[row['location']] += '\n>{i}_{g}\n{s}'.format(
            i=index, g=row['gene_name'], s=row['seq']
        )
        stats.setdefault(row['gene_name'], {'locations': []})
        stats[row['gene_name']]['locations'].append(row['location'])
    fname = {
        """3'UTR""": '3prime_utr.fa',
        """5'UTR""": '5prime_utr.fa',
        """CDS""": 'cds.fa',
        """ncRNA""": 'ncRNA.fa'
    }
    for label in seqs:
        with open('data/fasta/{f}'.format(f=fname[label]), 'w') as fh:
            fh.write(seqs[label])
    with open('data/fasta/all_locations_combined_fbfs.fa', 'w') as fh:
        fh.write("".join(seqs.values()))
    five_and_three = set()
    cds_and_three = set()
    five_cds_and_three = set()
    five_and_cds = set()
    only_three = set()
    only_cds = set()
    only_five = set()
    has_secondary = set()
    has_five = set()
    has_three = set()
    has_cds = set()
    for gene in stats:
        if ("""3'UTR""" in stats[gene]['locations']) and (
            """5'UTR""" in stats[gene]['locations']):
            five_and_three.add(gene)
        if ("""3'UTR""" in stats[gene]['locations']) and (
            """CDS""" in stats[gene]['locations']):
            cds_and_three.add(gene)
        if ("""3'UTR""" in stats[gene]['locations']) and (
            """5'UTR""" in stats[gene]['locations']) and (
            """CDS""" in stats[gene]['locations']
        ):
            five_cds_and_three.add(gene)
        if ("""5'UTR""" in stats[gene]['locations']) and (
            """CDS""" in stats[gene]['locations']):
            five_and_cds.add(gene)
        if "5'UTR" in stats[gene]['locations']:
            has_five.add(gene)
        if "3'UTR" in stats[gene]['locations']:
            has_three.add(gene)
        if "CDS" in stats[gene]['locations']:
            has_cds.add(gene)

        if set(stats[gene]['locations']) == set(["3'UTR"]):
            only_three.add(gene)
        if set(stats[gene]['locations']) == set(["5'UTR"]):
            only_five.add(gene)
        if set(stats[gene]['locations']) == set(["CDS"]):
            only_cds.add(gene)
        if len(stats[gene]['locations']) > 1:
            has_secondary.add(gene)
    print """
    Number of genes: {total}
    Number of genes with a 3' peak: {a_three}
    Number of genes with a CDS peak: {a_cds}
    Number of genes with a 5' peak: {a_five}
    Number of genes with only 3' peaks: {o_3}
    Number of genes with only CDS peaks: {o_cds}
    Number of genes with only 5' peaks: {o_5}
    Number of genes with 3' and 5' peaks: {a}
    Number of genes with 3' and CDS peaks: {b}
    Number of genes with 5' and CDS peaks: {five_and_cds}
    Number of genes with 5' and 3' and CDS peaks: {c}
    Number of genes with a secondary peak: {has_second}
    """.format(
        total=len(list(set(stats.keys()))),
        a_three=len(has_three),
        a_cds=len(has_cds),
        a_five=len(has_five),
        o_3=len(list(only_three)),
        o_cds=len(list(only_cds)),
        o_5=len(list(only_five)),
        a=len(list(five_and_three)),
        b=len(list(cds_and_three)),
        five_and_cds=len(five_and_cds),
        c=len(list(five_cds_and_three)),
        has_second=len(list(has_secondary))
    )
    labels = ['ncRNA', """3'UTR""", """5'UTR""", 'CDS']


def pick_ga(fname, ga):
    bfname = os.path.basename(fname)
    this_ga = None
    if bfname == 'combined_fbf1.txt':
        this_ga = ga['fbf1']
    if bfname == 'combined_fbf2.txt':
        this_ga = ga['fbf2']
    if bfname == 'combined_fbf.txt':
        this_ga = ga['both']
    return this_ga


def add_info_for_text(ga, args):
    for n in range(0, 14): print "~" * abs(7 - abs(n-7))
    print 'Adding info for input directories %s and %s.' % (
        args.five_reps, args.separate)
    for fname in glob.glob(args.five_reps + '/combined_fbf*.txt'):
        print "*" * 14 + '\n' + fname
        this_ga = pick_ga(fname, ga)
        if this_ga is None: continue
        peaks = pandas.read_csv(fname, sep='\t')
        for_every_row_get_heights_over_motifs(peaks, this_ga, ['tgt\w\w\wat', 'tgt[ag]aat'])
        get_info_from_combined_list_for_text(peaks, this_ga, fname, args.five_reps)

    get_info_from_separate_list_for_text(sep_dirname=args.separate, ga=ga)


def get_info_from_combined_list_for_text(peaks, this_ga, fname, outdir):
        """Write info for the main text.
        Write:
        1. Number of peaks and gene targets.
        2. Missing positives.
        3. How many positives are below rank 700.
        4. Fraction of peaks with or without an FBE.
        5. Number of genes with peaks in CDS vs 3'UTR, ect.
        6. Read coverage over 7 or 8mer motifs and t-test for the difference.
        """
        write_fasta_by_location(peaks)
        col_order = peaks.columns
        p = peaks.to_dict('records')
        for row in p:
            h = get_height(row, this_ga)
            if h is not None:
                row['height'] = h
        peaks = pandas.DataFrame(p)
        peaks.sort(columns='height', inplace=True, ascending=False)
        peaks = sort_cols(peaks)
        peaks['Rank'] = range(0, len(peaks.index))
        score_fbe(peaks)
        print score_positives(peaks)
        print "Number of peaks {p}. Number of gene targets {n}.".format(
            p=len(peaks.index),
            n=len(list(set(peaks['gene_name'].tolist())))
        )
        peaks.to_csv(outdir + '/%s' % os.path.basename(fname),  sep='\t', index=False)


def get_info_from_separate_list_for_text(sep_dirname=None, ga=None):
    """
    1. Score positives.
    2. Score FBE
    3. Write RIP shared/novel target numbers
    4. Write fastas of do/do not have FBE
    """
    peaks = {}
    try:
        peaks['fbf1'] = pandas.read_csv(sep_dirname + '/combined_fbf1.txt', sep='\t')
        peaks['fbf2'] = pandas.read_csv(sep_dirname + '/combined_fbf2.txt', sep='\t')
    except FileNotFoundError:
        print ('Was not able to read peaks files. Guessed {i} and {i2}.'.format(
            i=str(sep_dirname + 'combined_fbf1.txt'), i2=str(sep_dirname + 'combined_fbf2.txt')
        ))
        return

    print 'FBF-1'
    print score_positives(peaks['fbf1'])
    score_fbe(peaks['fbf1'])
    print 'FBF-2'
    print score_positives(peaks['fbf2'])
    score_fbe(peaks['fbf2'])
    num_peak_and_gene_numbers(peaks)
    num_shared = num_shared_targets(peaks)
    num_novel_vs_rip = find_num_novel_vs_rip(peaks)
    write_peaks_split_by_fbe(peaks['fbf1'], 'fbf1', sep_dirname)
    write_peaks_split_by_fbe(peaks['fbf2'], 'fbf2', sep_dirname)


def write_peaks_split_by_fbe(peaks, which_fbf, sep_dirname):
    peaks = peaks[peaks['biotype']=='protein_coding']
    with_fbe = peaks[peaks['has_fbe']==1]
    without_fbe = peaks[peaks['has_fbe']==0]
    if not os.path.exists('data/fasta/{f}/'.format(f=sep_dirname)):
        os.system('mkdir data/fasta/{f}/'.format(f=sep_dirname))
    li = ''
    for index, row in with_fbe.iterrows():
        li += '>{i}_{n}\n'.format(i=row['Rank'], n=row['gene_name'])
        li += row['seq'] + '\n'
    with open('data/fasta/{f}/{fbf}_{b}'.format(
            f=sep_dirname, fbf=which_fbf, b='with_fbe.fa'), 'w') as fh:
        fh.write(li)
    li = ''
    for index, row in without_fbe.iterrows():
        li += '>{i}_{n}\n'.format(i=row['Rank'], n=row['gene_name'])
        li += row['seq']+ '\n'
    with open('data/fasta/{f}/{fbf}_{b}'.format(
            f=sep_dirname, fbf=which_fbf, b='without_fbe.fa'), 'w') as fh:
            fh.write(li)
    li = ''
    for index, row in peaks.iterrows():
        li += '>{i}_{n}\n'.format(i=row['Rank'], n=row['gene_name'])
        li += row['seq']+ '\n'
    with open('data/fasta/{f}/{fbf}_{b}'.format(
            f=sep_dirname, fbf=which_fbf, b='all_peaks.fa'), 'w') as fh:
            fh.write(li)


def num_peak_and_gene_numbers(peaks):
    genes_1 = len(list(set(peaks['fbf1']['gene_name'])))
    genes_2 = len(list(set(peaks['fbf2']['gene_name'])))
    peaks_1 = len(peaks['fbf1'].index)
    peaks_2 = len(peaks['fbf2'].index)
    print """
    FBF-1 peaks: {p1} FBF-1 gene targets: {g1}
    FBF-2 peaks: {p2} FBF-2 gene targets: {g2}""".format(
        p1=peaks_1, p2=peaks_2, g1=genes_1, g2=genes_2
    )


def num_shared_targets(peaks):
    shared = set(peaks['fbf1']['gene_name'].tolist()) & set(
        peaks['fbf2']['gene_name'].tolist()
    )
    return len(list(shared))


def find_num_novel_vs_rip(peaks):
    rip, sam_ranks = compare_with_ripchip.get_ripchip_targets(get_sam_rank=True)
    either = set(peaks['fbf1']['gene_name'].tolist()) | set(
        peaks['fbf2']['gene_name'].tolist()
    )
    num_clip = len(list(either))
    num_rip = len(list(set(rip)))
    num_only_clip =len(list(
        either - set(rip)
    ))
    num_only_rip = len(list(
        set(rip) - either
    ))
    num_both = len(list(
        set(rip) & either
    ))
    num_either_rip_or_clip = len(list(
        set(rip) | either
    ))
    num_neither = 20000 - num_either_rip_or_clip
    fraction_of_clip = float(num_both)/float(num_clip)
    fraction_of_rip = float(num_both)/float(num_rip)
    print """Total clip: {tc}
    Total RIP: {tr}
    
    Only CLIP {oc}
    Only RIP: {orip}
    Both {b}
          """.format(
        tc=num_clip, tr=num_rip, oc=num_only_clip, orip=num_only_rip,
        b=num_both)
    print "New targets: %i" % num_only_clip
    print "Fraction of CLIP in RIP: %f" % float(fraction_of_clip)
    print "Fraction of RIP in CLIP: %f" % float(fraction_of_rip)

    res = scipy.stats.fisher_exact([
        [num_only_clip, num_both],
        [num_neither, num_only_rip]])
    print res
