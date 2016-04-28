"""
Create combined_fbf*.txt files that include both fbf1 and fbf2 read
information for each file of peaks.

1. Read in bedgraph information
2. Read in peaks lists
3. Apply read information to peaks
4. Output to unfiltered/
5. Filter and output to filtered/

# config.ini format:
[library]
top: /groups/Kimble/Aman Prasad/clip2/
bedgraphs_folder: %(top)s/bedgraph_unnorm/
read_beds: %(top)s/bed_collapsed/
fbf1_peaks: %(top)s/fbf1/peaks/combined_fbf1/null_hyp_4.txt
fbf2_peaks: %(top)s/fbf2/peaks/combined_fbf2/null_hyp_4.txt
fbf1_rep_1: %(top)s/fbf1/bedgraph_norm/exp_fbf1_TGGC
fbf1_rep_2: %(top)s/fbf1/bedgraph_norm/exp_fbf1_GGTT
fbf1_rep_3: %(top)s/fbf1/bedgraph_norm/exp_fbf1_CGGA
fbf2_rep_1: %(top)s/fbf2/bedgraph_norm/exp_fbf2_TGGC
fbf2_rep_2: %(top)s/fbf2/bedgraph_norm/exp_fbf2_GGTT
fbf2_rep_3: %(top)s/fbf2/bedgraph_norm/exp_fbf2_CGGA
fbf1_n2: %(top)s/fbf1/bedgraph_norm/control_n2
fbf2_n2: %(top)s/fbf2/bedgraph_norm/control_n2
rna_seq: %(top)s/bedgraph_norm/rna_seq
"""
import HTSeq
import re
import sys
import os
import pandas
import argparse
import config
import numpy as np

from peaks import peaks
        

def load_bedgraph(fname):
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    plus_file = fname.partition('.wig')[0] + '_+.wig'
    add_strand_to_ga_from_bedgraph_file(plus_file, ga, '+')
    minus_file = fname.partition('.wig')[0] + '_-.wig'
    add_strand_to_ga_from_bedgraph_file(minus_file, ga, '-')
    return ga


def add_strand_to_ga_from_bedgraph_file(fname, ga, strand):
    with open(fname, 'r') as f:
        next(f)
        for line in f:
            #if line[0] != 'M': continue
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), strand)] = float(s[3])
    return ga


def get_val(ga, iv):
    return np.max(np.fromiter(ga[iv], dtype=np.float))


def get_pos_of_max_coverage(ga, iv):
    index = np.argmax(np.fromiter(ga[iv], dtype=np.float))
    return iv.start + index


def add_heights_to_peak_file(peak_fname, bed_ga_dict):
    try:
        peaks = pandas.read_csv(peak_fname, sep='\t')
    except:
        print "Could not open %s" % peak_fname
        return


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
    open(filename, 'w').write(separate)


def do_we_skip_adding_reads(fbf1, fbf2, args, reps):
    skip = False
    if (fbf1.has_columns(reps) and fbf2.has_columns(reps)):
        if args.no_ui:
            return False
        answ = raw_input('''Both peaks files appear to have reads already. \
Use those? [y/n]?''').lower()
        if answ[0] == 'y':
            skip = True
        elif answ[0] == 'n':
            skip = False
        else:
            print  "Invalid response."
            sys.exit()
    return skip


def do_we_skip_adding_reads_to_replicates(args):
    if args.no_ui: return False
    if raw_input(
        'Add info to individual replicates? [y/n]:').lower() == 'y':
        return False
    else: return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
''')
    parser.add_argument('-c', '--config_ini')
    parser.add_argument('-w', '--write_example_config',
                        default=False, action='store_true')
    parser.add_argument('-r', '--ratio_cutoff', default='5',
                        help='Enrichment ratio cutoff.')
    parser.add_argument('-o', '--output', default='output.txt')
    parser.add_argument('-v', '--no_ui', action='store_true',
                        default=False)
    args = parser.parse_args()
    if args.write_example_config:
        outfname = 'tmp_example_config.ini'
        if os.path.exists(outfname):
            print "Could not write example config.ini because \
{a} file exists.".format(a=outfname)
            sys.exit()
        else: write_config(outfname)
        sys.exit()
    args.ratio_cutoff = float(args.ratio_cutoff)
    lib = config.config(filepath=args.config_ini)
    fbf1 = peaks(file=lib['fbf1_peaks'], name='fbf1')
    fbf2 = peaks(file=lib['fbf2_peaks'], name='fbf2')
    reps_labels = ['fbf1_rep_1', 'fbf1_rep_2', 'fbf1_rep_3',
        'fbf2_rep_1', 'fbf2_rep_2', 'fbf2_rep_3',
            'rna_seq', 'rna_seq_modencode', 'rna_seq_oo', 'fbf1_n2', 'fbf2_n2']
    rep = {}
    for rootdir in [lib['fbf1_peak_reps'], lib['fbf2_peak_reps']]:
        for subdir, dirs, files in os.walk(rootdir):
            for file_basename in files:
                file_path = os.path.join(subdir, file_basename)
                rep[file_path] = peaks(file=file_path, name=file_basename)
                print rep[file_path]
    skip_combined = do_we_skip_adding_reads(fbf1, fbf2, args, reps_labels)
    skip_reps = do_we_skip_adding_reads_to_replicates(args)
    if not skip_combined:
        ga = {}
        for (name, fname) in[(x, lib[x]) for x in reps_labels]:
            ga[name] = load_bedgraph(fname)
            fbf1.add_reads(ga=ga[name], name=name)
            fbf2.add_reads(ga=ga[name], name=name)
    if not skip_reps:
        for rootdir in [lib['fbf1_peak_reps'], lib['fbf2_peak_reps']]:
            for subdir, dirs, files in os.walk(rootdir):
                for file_basename in files:
                    file_path = os.path.join(subdir, file_basename)
                    for (name, fname) in[(x, lib[x]) for x in reps_labels]:
                        rep[file_path].add_reads(ga=ga[name], name=name)
                    print "Start:"
                    print rep[file_path]
                    rep[file_path].set_sum(
                        to_sum=['fbf1_rep_1', 'fbf1_rep_2', 'fbf1_rep_3'],
                        summed_col='fbf1_reads')
                    rep[file_path].set_sum(
                        to_sum=['fbf2_rep_1', 'fbf2_rep_2', 'fbf2_rep_3'],
                        summed_col='fbf2_reads')
                    if re.search('fbf1', name):
                        rep[file_path].set_ratio(col1='fbf1_reads', col2='fbf1_n2')
                    if re.search('fbf2', name):
                        rep[file_path].set_ratio(col1='fbf2_reads', col2='fbf2_n2')
                    rep[file_path].write_table(file_path)
                    print "End:"
                    print rep[file_path]
    fbf1.set_sum(
        to_sum=['fbf1_rep_1', 'fbf1_rep_2', 'fbf1_rep_3'],
        summed_col='fbf1_reads')
    fbf1.set_sum(
        to_sum=['fbf2_rep_1', 'fbf2_rep_2', 'fbf2_rep_3'],
        summed_col='fbf2_reads')
    fbf2.set_sum(
        to_sum=['fbf1_rep_1', 'fbf1_rep_2', 'fbf1_rep_3'],
        summed_col='fbf1_reads')
    fbf2.set_sum(
        to_sum=['fbf2_rep_1', 'fbf2_rep_2', 'fbf2_rep_3'],
        summed_col='fbf2_reads')
    fbf1.set_ratio(col1='fbf1_reads', col2='fbf1_n2')
    fbf2.set_ratio(col1='fbf2_reads', col2='fbf2_n2')
    filt1 = fbf1.get_filtered_obj(col='ratio', cutoff=args.ratio_cutoff)
    filt2 = fbf2.get_filtered_obj(col='ratio', cutoff=args.ratio_cutoff)
    fbf1.write_table('unfiltered_separate/combined_fbf1.txt')
    fbf2.write_table('unfiltered_separate/combined_fbf2.txt')
    filt1.write_table('filtered_separate/combined_fbf1.txt')
    filt2.write_table('filtered_separate/combined_fbf2.txt')
