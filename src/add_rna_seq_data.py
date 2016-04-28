import HTSeq
import re
import sys
import os
import pandas
import argparse
import config
import numpy as np
sys.path.insert(
    0, '/groups/Kimble/Aman Prasad/clip2/src/')
from peaks import peaks
from combine_peak_information import load_bedgraph
del sys.path[0]


def write_config(filename):
    # config.ini format:
    separate ='''[library]
top: /groups/Kimble/Aman Prasad/clip2/
bedgraphs_folder: %(top)s/bedgraph_unnorm/
read_beds: %(top)s/bed_collapsed/
fbf1_peak_reps: %(top)s/fbf1/peaks/
fbf2_peak_reps: %(top)s/fbf2/peaks/
fbf1_peaks: %(top)s/fbf1/peaks/combined_fbf1/null_hyp_4.txt
fbf2_peaks: %(top)s/fbf2/peaks/combined_fbf2/null_hyp_4.txt
# Bedgraph files.
fbf1_rep_1: %(top)s/fbf1/bedgraph_norm/exp_fbf1_TGGC
fbf1_rep_2: %(top)s/fbf1/bedgraph_norm/exp_fbf1_GGTT
fbf1_rep_3: %(top)s/fbf1/bedgraph_norm/exp_fbf1_CGGA
fbf2_rep_1: %(top)s/fbf2/bedgraph_norm/exp_fbf2_TGGC
fbf2_rep_2: %(top)s/fbf2/bedgraph_norm/exp_fbf2_GGTT
fbf2_rep_3: %(top)s/fbf2/bedgraph_norm/exp_fbf2_CGGA
fbf1_n2: %(top)s/fbf1/bedgraph_norm/control_n2
fbf2_n2: %(top)s/fbf2/bedgraph_norm/control_n2
rna_seq: %(top)s/bedgraph_norm/rna_seq
rna_seq_modencode: %(top)s/bedgraph_norm/rna_seq_mod4594
rna_seq_oo: %(top)s/bedgraph_norm/rna_seq_oo_srr1263137

# Bed files.
fbf1_rep_1_bed: %(top)s/fbf1/bed_collapsed/exp_fbf1_TGGC
fbf1_rep_2_bed: %(top)s/fbf1/bed_collapsed/exp_fbf1_GGTT
fbf1_rep_3_bed: %(top)s/fbf1/bed_collapsed/exp_fbf1_CGGA
fbf2_rep_1_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_TGGC
fbf2_rep_2_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_GGTT
fbf2_rep_3_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_CGGA
fbf1_n2_bed: %(top)s/bed_collapsed/fbf1_n2
fbf2_n2_bed: %(top)s/bed_collapsed/fbf2_n2
rna_seq_bed: %(top)s/bed_collapsed/rna_seq
rna_seq_modencode_bed: %(top)s/bed_collapsed/rna_seq_mod4594
rna_seq_oo_bed: %(top)s/bed_collapsed/rna_seq_oo_srr1263137
'''
    same_n2 ='''[library]
top: /groups/Kimble/Aman Prasad/clip2/
bedgraphs_folder: %(top)s/bedgraph_unnorm/
read_beds: %(top)s/bed_collapsed/
fbf1_peak_reps: %(top)s/fbf1_to_fbf2_n2/peaks/
fbf2_peak_reps: %(top)s/fbf2/peaks/
fbf1_peaks: %(top)s/fbf1_to_fbf2_n2/peaks/combined_fbf1_to_fbf2_n2/null_hyp_4.txt
fbf2_peaks: %(top)s/fbf2/peaks/combined_fbf2/null_hyp_4.txt
# Bedgraph files.
fbf1_rep_1: %(top)s/fbf1_to_fbf2_n2/bedgraph_norm/exp_fbf1_TGGC
fbf1_rep_2: %(top)s/fbf1_to_fbf2_n2/bedgraph_norm/exp_fbf1_GGTT
fbf1_rep_3: %(top)s/fbf1_to_fbf2_n2/bedgraph_norm/exp_fbf1_CGGA
fbf2_rep_1: %(top)s/fbf2/bedgraph_norm/exp_fbf2_TGGC
fbf2_rep_2: %(top)s/fbf2/bedgraph_norm/exp_fbf2_GGTT
fbf2_rep_3: %(top)s/fbf2/bedgraph_norm/exp_fbf2_CGGA
fbf1_n2: %(top)s/fbf1_to_fbf2_n2/bedgraph_norm/control_n2
fbf2_n2: %(top)s/fbf2/bedgraph_norm/control_n2
rna_seq: %(top)s/bedgraph_norm/rna_seq
rna_seq_modencode: %(top)s/bedgraph_norm/rna_seq_mod4594
rna_seq_oo: %(top)s/bedgraph_norm/rna_seq_oo_srr1263137
# Bed files.
fbf1_rep_1_bed: %(top)s/fbf1_to_fbf2_n2/bed_collapsed/exp_fbf1_TGGC
fbf1_rep_2_bed: %(top)s/fbf1_to_fbf2_n2/bed_collapsed/exp_fbf1_GGTT
fbf1_rep_3_bed: %(top)s/fbf1_to_fbf2_n2/bed_collapsed/exp_fbf1_CGGA
fbf2_rep_1_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_TGGC
fbf2_rep_2_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_GGTT
fbf2_rep_3_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_CGGA
fbf1_n2_bed: %(top)s/bed_collapsed/fbf1_n2
fbf2_n2_bed: %(top)s/bed_collapsed/fbf2_n2
rna_seq_bed: %(top)s/bed_collapsed/rna_seq
rna_seq_modencode_bed: %(top)s/bed_collapsed/rna_seq_mod4594
rna_seq_oo_bed: %(top)s/bed_collapsed/rna_seq_oo_srr1263137
'''
    open(filename, 'w').write(same_n2)


def to_rpkm(x, _dict):
    if x not in _dict:
        return ''
    else:
        return float(_dict[x])

def add_abundance(_peak):
    print _peak.data.gene_name
    print '*'
    print _peak
    df = pandas.read_csv('lib/ortiz/DESeq_genes_in_gonad.txt', sep='\t')
    df['rpkm'] = df['Expression in  fog-2 gonads (RPKM)']
    gene_to_rpkm = dict(zip(df['Gene name'].tolist(), df.rpkm))
    _peak.data['RNA abundance'] = [to_rpkm(x, gene_to_rpkm) \
                                   for x in _peak.data.gene_name]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
''')
    parser.add_argument('-c', '--config_ini')
    parser.add_argument('-w', '--write_example_config',
                        default=False, action='store_true')
    parser.add_argument('-o', '--output', default='output.txt')
    args = parser.parse_args()
    if args.write_example_config:
        outfname = 'tmp_example_config.ini'
        if os.path.exists(outfname):
            print "Could not write example config.ini because \
{a} file exists.".format(a=outfname)
            sys.exit()
        else: write_config(outfname)
        sys.exit()
    lib = config.config(filepath=args.config_ini)
    lib['fbf1_filt'] = 'filtered/combined_fbf1.txt'
    lib['fbf1_separate_filt'] = 'filtered_separate/combined_fbf1.txt'
    lib['fbf1_separate_unfilt'] = 'unfiltered_separate/combined_fbf1.txt'
    lib['fbf2_filt'] = 'filtered/combined_fbf2.txt'
    lib['fbf1_unfilt'] = 'unfiltered/combined_fbf1.txt'
    lib['fbf2_unfilt'] = 'unfiltered/combined_fbf2.txt'
    lib['both_filt'] = 'filtered_five_reps/combined_fbf.txt'
    rep = {}
    ga = {}
    for rootdir in [lib['fbf1_peak_reps'], lib['fbf2_peak_reps']]:
        for subdir, dirs, files in os.walk(rootdir):
            for file_basename in files:
                file_path = os.path.join(subdir, file_basename)
                rep[file_path] = peaks(file=file_path, name=file_basename)
                print rep[file_path]
    fbf1 = peaks(file=lib['fbf1_peaks'], name='fbf1')
    fbf2 = peaks(file=lib['fbf2_peaks'], name='fbf2')
    both = peaks(file='filtered_five_reps/combined_fbf.txt', name='5reps')
    fbf1_filt = peaks(file=lib['fbf1_filt'], name='fbf1_filt')
    fbf1_separate_filt = peaks(file=lib['fbf1_separate_filt'], name='fbf1_separate_filt')
    fbf1_separate_unfilt = peaks(file=lib['fbf1_separate_unfilt'], name='fbf1_separate_unfilt')
    fbf2_filt = peaks(file=lib['fbf2_filt'], name='fbf2_filt')
    fbf1_unfilt = peaks(file=lib['fbf1_unfilt'], name='fbf1_unfilt')
    fbf2_unfilt = peaks(file=lib['fbf2_unfilt'], name='fbf2_unfilt')
    both_filt = peaks(file=lib['both_filt'], name='both_filt')
    for x in [fbf1, fbf2, both, fbf1_filt, fbf2_filt, fbf1_unfilt,
              fbf2_unfilt, both_filt,
              fbf1_separate_filt, fbf1_separate_unfilt]:
        add_abundance(x)
    for x in rep: add_abundance(rep[x])
    reps_labels = ['rna_seq_modencode', 'rna_seq_oo']
    for (name, fname) in[(x, lib[x]) for x in reps_labels]:
        ga[name] = load_bedgraph(fname)
        for x in [fbf1, fbf2, both, fbf1_filt, fbf2_filt, fbf1_unfilt,
              fbf2_unfilt, both_filt, fbf1_separate_filt, fbf1_separate_unfilt]:
            x.add_reads(ga=ga[name], name=name)
        for f in rep:
            rep[f].add_reads(ga=ga[name], name=name)
    for x in [fbf1, fbf2, both, fbf1_filt, fbf2_filt, fbf1_unfilt,
              fbf2_unfilt, both_filt, fbf1_separate_filt, fbf1_separate_unfilt]:
        x.write_table(x.file)
    for x in rep:
        rep[x].write_table(rep[x].file)
