import os
import glob
import pandas
import _filter
import config
import sys
import argparse
import re

def paths(config):
    if config is not None:
        non_control = {}
        bedfiles = {'control': "{a}/{b}.wig".format(
            a=os.path.dirname(config['clip_replicate'][0]),
            b=os.path.basename(config['neg_ip_filename']).rstrip('.bed'))}
        for x in config['clip_replicate']:
            name = os.path.basename(x).rstrip('.bed').rstrip('.wig')
            bedfiles[name] = x
            non_control[name] = x
    print bedfiles
    return bedfiles, non_control

def mk(top):
    if not os.path.exists(top): os.system('mkdir ' + top)

def fill_in_params(topdir):
    reads_label = None
    m = re.search('(fbf\d)', topdir)
    look_for = 'fbf1'
    if m is not None:
        answ = raw_input(
'''This looks like {a}. Label the per million reads averaged \
across replicates as {a}_reads? (Y/N):'''.format(a=m.group(1)))
        if answ.upper() == 'Y':
            reads_label = '{a}.format(a=m.group(1))
            print "will fill in as {a}_reads" + reads_label
        if m.group(1) == 'fbf1':
            look_for = 'fbf2'
    guess_path = '/groups/Kimble/Aman Prasad/clip2/unfiltered/' + \
                      'combined_' + look_for + '.txt'
    answ = raw_input(
'''Are we filling in 1. replicates or 2, both FBF-1 and FBF-2 (+N2s)?
Replicates: {a}
Other FBF: {b}
We will use the exp and control columns in the case of the later.
[1/2]>
'''.format(a=bedfiles, b=guess_path))
    if answ == '1':
        filling_in = 'replicates'
    else:
        filling_in = 'other fbf'
        if os.path.exists(guess_path):
            try:
                test = pandas.read_csv(guess_path, sep='\t')
                print "(Testing {a}...{b} lines OK.)".format(
                    a=guess_path, b=[len(test.exp), len(test.control)])
            except:
                print "Problem with {z}. Aborting...".format(
                    z=guess_path)
                sys.exit()
    if reads_label == 'fbf1':
        other_label = 'fbf2'
    if reads_label == 'fbf2':
        other_label = 'fbf1'
    return reads_label, filling_in, other_label

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config')
args = parser.parse_args()

def fill_in(peaks, reads_label, other_file, other_label):
    if reads_label is not None:
        peaks[reads_label + '_reads'] = peaks['exp']
        peaks['n2_for_' + reads_label + '_reads'] = peaks['control']
        peaks[reads_label + '_n2'] = peaks['control']
    other = pandas.read_csv(guess_path, sep='\t')

lib = config.config(args.config)
if 'experiment_name' in lib:
    topdir = lib['experiment_name'] + '/peaks/' # 'with_nums/'
else:
    topdir = 'fbf1/peaks/'
if not os.path.exists(topdir):
    print "Could not find %s." % topdir

bedfiles, non_control = paths(lib)
reads_label, other_file, other_label = fill_in_params(topdir)
ga = {}
print "Loading bedgraphs..."
for name in bedfiles:
    print "Loading %s" % name
    ga[name] = _filter.load_bedgraph_file(bedfiles[name])
for fname in glob.glob(topdir + '/*/*txt'):
    print "Adding heights to %s" % fname
    peaks = _filter.add_heights_to_peak_file(fname, ga)
    if len(peaks.index) == 0:
        print "%s empty..." % fname
        continue
    peaks = _filter.normalize_height_columns(bedfiles, lib, peaks)
    peaks = _filter.add_sum_and_ratio_columns(peaks)
    keep_col = set([x for x in peaks.columns if not (
        re.search('local', x) or re.search('_norm', x) \
        or re.search('poisson', x) or re.search('pvalues', x) or \
        re.search('exons', x))])
    peaks = peaks[list(keep_col)]
    fill_in(peaks, reads_label, other_file, other_label)
    top = 'with_nums/'
    mk(top)
    subdir = os.path.basename(os.path.dirname(os.path.realpath(fname)))
    mk(top + subdir)
    outfname = top + '/' + subdir + '/' + os.path.basename(fname)
#    if not os.path.exists(outfname):
    print "to " + outfname
    peaks.to_csv(outfname, sep='\t', index=False)
