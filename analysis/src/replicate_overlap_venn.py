import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib
import sys
from matplotlib_venn import *
import os
import collections


def make_venn_of_combined(peaks, label='', subdir=''):
    targets = {}
    in_any = set()

    for replicate in peaks:
        targets[replicate] = set(peaks[replicate]['gene_name'])
        in_any = in_any | targets[replicate]
        print "{r}: {n} targets.".format(r=replicate, n=len(targets[replicate]))
    in_all = in_any
    for replicate in peaks:
        in_all = in_all & targets[replicate]
    print "In any replicate {N}. In all replicates {a}".format(
        N=len(in_any), a=len(in_all)
    )
    n_unique_to = collections.defaultdict()
    for replicate in peaks:
        n_unique_to[replicate] = len(targets[replicate] - in_all)
        print "Unique to {r}: {n}. {n} + {a} = {na}".format(
            r=replicate, n=n_unique_to[replicate], a=len(in_all),
            na=len(in_all) + n_unique_to[replicate]
        )
    plt.clf()
    overlaps = venn2(
        [x for x in targets.values()],
        (os.path.basename(x).partition('.txt')[0] for x in targets))
    if not os.path.exists(subdir): os.system('mkdir %s' % subdir)
    plt.savefig('%s/fbf1_vs_fbf2_overlap%s.pdf' % (subdir, label), format='pdf')
    plt.clf()

def _mk(dirname):
    if not os.path.exists(dirname):
        os.system('mkdir ' + dirname)

def write_tables_of_joint_and_unique(peaks, subdir=''):
    targets = {}
    subdir = re.sub('figs', '', subdir)
    for exp in peaks:
        exp_base = os.path.basename(exp).partition('.txt')[0]
        targets[exp_base] = set(peaks[exp]['gene_name'])
    out_dir = 'tables/' + subdir
    _mk('tables')
    _mk(out_dir)
    a_key = targets.keys()[0]
    b_key = targets.keys()[1]
    a_only_filename = out_dir + '/%s_only_targets.txt' % a_key
    b_only_filename = out_dir + '/%s_only_targets.txt' % b_key
    both_filename = out_dir + '/joint_targets.txt'
    print "\tWriting tables to %s/..." % out_dir
    with open(a_only_filename, 'w') as f:
        a_only_targets = targets[a_key] - targets[b_key]
        f.write(write_as_line(a_only_targets))
    with open(b_only_filename, 'w') as f:
        b_only_targets = targets[b_key] - targets[a_key]
        f.write(write_as_line(b_only_targets))
    with open(both_filename, 'w') as f:
        joint_targets = targets[b_key] & targets[a_key]
        f.write(write_as_line(joint_targets))              


def write_as_line(set_of_targs):
    li = ''
    for gene in set_of_targs:
        li += '%s\n' % gene
    return li        


def make_venn_of_replicates(peaks, label='no_label', subdir=''):
    targets = {}
    in_any = set()
    for replicate in peaks:
        targets[replicate] = set(peaks[replicate]['gene_name'].dropna())
        in_any = in_any | targets[replicate]
    #print targets
    in_all = in_any
    for replicate in peaks:
        in_all = in_all & targets[replicate]
    plt.clf()
    #for t in targets: targets[t] = list(targets[t])
    labels = [os.path.basename(x).partition('.txt')[0] for x in targets]
    labels = tuple(labels)
    overlaps = venn3(
        [targets[x] for x in targets],
        (x for x in targets))
    if not os.path.exists(subdir): os.system('mkdir %s' % subdir)
    print label
    print ">=1/3 reps have {i} targets/".format(i=len(in_any))
    print "3/3 reps have {i} targets/".format(i=len(in_all))
    print "\tsaving to ",
    print '%s/replicate_overlap_%s.pdf' % (subdir, label)
    plt.savefig('%s/replicate_overlap_%s.pdf' % (subdir, label), format='pdf')
    plt.clf()


def huh(peaks, label=''):
    targets = {}
    in_any = set()
    for replicate in peaks:
        targets[replicate] = set(peaks[replicate]['gene_name'].dropna())
        in_any = in_any | targets[replicate]
    #print targets
    in_all = in_any
    for replicate in peaks:
        in_all = in_all & targets[replicate]
    import pandas
    if label == 'fbf1':
        df = pandas.read_csv('filtered/combined_fbf1.txt', sep='\t')
    if label == 'fbf2':
        df = pandas.read_csv('filtered/combined_fbf2.txt', sep='\t')
    expected = set(df.gene_name)
    not_expected =  expected - in_all
    print "In target list but not 3/3 Venn: %i" % len(not_expected)
    print not_expected
    print "--" * 10

def survive_filter(peaks):
    targets = {}
    in_any = set()
    for replicate in peaks:
        targets[replicate] = set(peaks[replicate]['gene_name'].dropna())
        in_any = in_any | targets[replicate]
    return in_any
    
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-r1', '--replicates_dir_fbf1', default=None,
                    help='A directory of biological replicates for fbf1')
parser.add_argument('-r2', '--replicates_dir_fbf2', default=None,
                    help='A directory of biological replicates for fbf2')
parser.add_argument('-c', '--combined_dir', default=None)
parser.add_argument('-o', '--output_dir', default='figs/')
args = parser.parse_args()
#gtf_sep_cols = pandas.read_csv('lib/gtf_with_names_column.txt', sep='\t')
#top_level_dir = 'pypeaks_fdr5_negip_local/'
if (args.replicates_dir_fbf1 is not None) and (
    args.replicates_dir_fbf2 is not None):
    fbf1 = {}
    filtered_fbf1 = {}
    fbf2 = {}
    filtered_fbf2 = {}
    for filename in glob.glob(args.replicates_dir_fbf1 + '/*fbf1*.txt'):
        fbf1[filename] = pandas.read_csv(filename, sep='\t')
    for filename in glob.glob(args.replicates_dir_fbf2 + '/*fbf2*.txt'):
        fbf2[filename] = pandas.read_csv(filename, sep='\t')
    #try:
    for filename in glob.glob('biological_reps_filtered/*fbf1*.txt'):
        filtered_fbf1[filename] = pandas.read_csv(filename, sep='\t')
    for filename in glob.glob('biological_reps_filtered/*fbf2*.txt'):
        print filename
        filtered_fbf2[filename] = pandas.read_csv(filename, sep='\t')
    print filtered_fbf1
    fbf1_survive_filter = survive_filter(filtered_fbf1)
    fbf2_survive_filter = survive_filter(filtered_fbf2)
    print "\n" * 20
    if 'top-1' in fbf2_survive_filter:
        print 'in'
    else: print 'out'
    #except:
    #    print "No filtered dir."
    print "\tMaking venn of replicates from %i files..." % (
        int(len(fbf1) + len(fbf2)))
    make_venn_of_replicates(fbf1, label='unfilt_fbf1',
                            subdir=args.output_dir)
    huh(fbf1, label='fbf1')
    make_venn_of_replicates(fbf2, label='unfilt_fbf2',
                            subdir=args.output_dir)
    huh(fbf2, label='fbf2')
    skip = '''
    for t in fbf1:
        fbf1[t]['survive_filter'] = [(x in fbf1_survive_filter) for x in fbf1[t].gene_name.dropna()]
        fbf1[t] = fbf1[t][fbf1[t]['survive_filter']]
    for t in fbf2:
        print "===FBF2===="
        print len(fbf2[t].index)
        fbf2[t]['survive_filter'] = [(x in fbf2_survive_filter) for x in fbf2[t].gene_name.dropna()]
        print fbf2[t]['survive_filter'].value_counts()
        fbf2[t] = fbf2[t][fbf2[t]['survive_filter']]
        print '>'
        print len(fbf2[t].index)'''
    make_venn_of_replicates(filtered_fbf1, label='fbf1',
                            subdir=args.output_dir)
    huh(fbf1, label='fbf1')
    make_venn_of_replicates(filtered_fbf2, label='fbf2',
                            subdir=args.output_dir)
    huh(fbf2, label='fbf2')
if args.combined_dir is not None:
    combined = {}
    for filename in glob.glob(args.combined_dir + '/combined*.txt'):
        combined[filename] = pandas.read_csv(filename, sep='\t')
    print "\tMaking venn of combined targets from %i files..." % len(combined)
    make_venn_of_combined(combined, subdir=args.output_dir)
    write_tables_of_joint_and_unique(combined,
                                     subdir=args.output_dir)
