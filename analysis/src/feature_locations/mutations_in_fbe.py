import re
import sys
import os
import random
import matplotlib.pyplot as plt
import numpy as np


def mutation_dict(motif='nntgtnnnatn'):
    mutations_regex = [
    ]
    bases = ['a', 't', 'c', 'g']
    mutations = {}
    for pos in range(len(motif)):
        seq = motif.lower()
        for base in bases:
            seq = list(seq)
            seq[pos] = base
            seq = ''.join(seq)
            tup = (seq, base, pos)
            mutations[tup] = seq.replace('n', r'''\w''')
            mutations[tup] = '.*' + mutations[tup] + '.*'
    return mutations


def read_fasta(fastq_filename):
    fasta_lines = []
    with open(fastq_filename, 'r') as f:
        for li in f:
            if re.match('\A>.*', li) is None:
                fasta_lines.append(li.rstrip('\n'))
    print fasta_lines[0]
    return fasta_lines


def scramble(fasta_lines):
    scrambled = []
    for line in fasta_lines:
        scrambled.append(''.join(random.sample(line, len(line))))
    return scrambled


def mutations(fasta_filename='data/fasta/subpeaks/highest_both.fa',
              output_name='figs/Bernstein.pdf',
              base_motif='NNUGUNNNAUN'):
    base_motif = base_motif.lower()
    base_motif = re.sub('u', 't', base_motif)
    print base_motif
    fasta_lines_pos = read_fasta(fasta_filename)
    mutations = mutation_dict(motif=base_motif)
    enrichment = dict()
    for tup in mutations:
        regex_pat = mutations[tup]
        enrichment[tup] = get_enrichment(regex_pat, fasta_lines_pos)
    make_fig(enrichment, base_motif=base_motif, output_name=output_name)


def make_fig(table, base_motif='NNUGUNNNAUN', output_name='figs/Bernstein.pdf'):
    assert isinstance(table, type({}))
    # table.keys() are tuples: (name, base, position in fbe)
    n_groups = len(list(
        set([x[2] for x in table.keys()])
    ))
    in_order = sorted(table.keys(), key=lambda tup: int(tup[2]))
    print in_order
    a_ratios = [table[x] for x in in_order if x[1]=='a']
    t_ratios = [table[x] for x in in_order if x[1]=='t']
    c_ratios = [table[x] for x in in_order if x[1]=='c']
    g_ratios = [table[x] for x in in_order if x[1]=='g']
    x_locs = np.arange(n_groups)
    width = 0.15
    fig, ax = plt.subplots()
    rects1 = ax.bar(x_locs, a_ratios, width, color='k', linewidth=0)
    rects2 = ax.bar(x_locs + width, t_ratios, width, color='r', linewidth=0)
    rects3 = ax.bar(x_locs + 2*width, c_ratios, width, color='b', linewidth=0)
    rects4 = ax.bar(x_locs + 3*width, g_ratios, width, color='g', linewidth=0)
    x_tick_locs = []
    x_tick_labels = []
    for loc in x_locs:
        x_tick_locs.append(loc)
        x_tick_labels.append('a')
        x_tick_locs.append(loc+width)
        x_tick_labels.append('t')
        x_tick_locs.append(loc+2*width)
        x_tick_labels.append('c')
        x_tick_locs.append(loc+3*width)
        x_tick_labels.append('g')
    plt.tick_params(
        axis='x',
        which='both',
        bottom='off',
        top='off',
        labelbottom='on'
    )
    ax.set_xlim([-0.05, x_locs[-1] + 5*width + 0.05])
    ax.set_xticks(x_tick_locs)
    ax.set_xticklabels(x_tick_labels)
    ax.set_xticks(x_locs + 1.5*width)
    ax.set_xticklabels([str(x)+'\n{z}'.format(z=list(base_motif)[x]) for x in range(n_groups)])
    legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),
              ('A', "U", "C", "G"),
              loc='upper center', ncol=4)
    ax.legend(bbox_to_anchor=(0.7, 1.1),
               #bbox_transform=plt.gcf().transFigure,
               linestyle='None',
               linewidths=0)
    #for label in ax.legend.get_texts():
    #    label.set_fontsize(8)
    frame = legend.get_frame()
    frame.set_linewidth(0)
    ax.set_ylabel('Fold enrichment in FBF peaks')
    ax.set_xlabel('Position of mutation')
    plt.savefig(output_name, format='pdf')


def get_enrichment(regex_pat, fasta_lines_pos):
    regex_obj = re.compile(regex_pat, re.IGNORECASE).match
    positive_seqs = score_in_fasta(fasta_lines_pos, regex_obj)
    negative_seqs = get_average_score_after_scrambling(
        fasta_lines_pos, regex_obj, num_permutations=1000)

    print "Positive sequences: with %s: %i, without: %i" % (
        regex_pat, positive_seqs[0], positive_seqs[1]
    )
    print "Negative sequences: with %s: %f, without: %f" % (
        regex_pat, negative_seqs[0], negative_seqs[1]
    )
    frac_pos = float(positive_seqs[0])/float(sum(positive_seqs))
    frac_neg = float(negative_seqs[0])/float(sum(negative_seqs))
    if frac_neg == 0:
        return -1
    return frac_pos/frac_neg


def score_in_fasta(fasta, regex_obj):
    score = [0, 0]  # with, without.
    for line in fasta:
        if regex_obj(line):
            score[0] += 1
        else:
            score[1] += 1
    return score


def get_average_score_after_scrambling(lines, regex_obj, num_permutations=1000):
    results = []
    for i in range(num_permutations):
        mutated_seqs = scramble(lines)
        results.append(np.array(
            score_in_fasta(mutated_seqs, regex_obj)
        ))
    results = np.array(results)
    # print results
    averages = np.average(results, axis=0)
    print averages
    return averages


if __name__ == '__main__':
    mutations(sys.argv[1], base_motif='NNUGUGNAUNN', output_name='figs/Bernstein_7mer.pdf')
    mutations(sys.argv[1], base_motif='NNUGUNNNAUN', output_name='figs/Bernstein_fbe.pdf')
