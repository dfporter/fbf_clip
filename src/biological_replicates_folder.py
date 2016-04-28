import glob
import os
import sys
sys.path.insert(
    0, './analysis/src/')
from peaks import peaks
import re

def scribe(cmd):
    os.system('echo')
    os.system('echo ' + cmd)
    os.system(cmd)


def _mk(name):
    if not os.path.exists(name):
        os.system('mkdir ' + name)

def run(top, cutoff=3):
    trans = {
        'exp_fbf1_TGGC': 'fbf1_rep_1',
        'exp_fbf1_GGTT': 'fbf1_rep_2',
        'exp_fbf1_CGGA': 'fbf1_rep_3',
        'exp_fbf2_TGGC': 'fbf2_rep_1',
        'exp_fbf2_GGTT': 'fbf2_rep_2',
        'exp_fbf2_CGGA': 'fbf2_rep_3',
        }

    for fname in glob.glob(top + '/peaks/exp_*/null_hyp_4.txt'):
        dname = os.path.basename(os.path.dirname(fname))
        print fname
        print trans[dname]
        p = peaks(file=fname)
        p.set_sum(
            to_sum=['fbf1_rep_1', 'fbf1_rep_2', 'fbf1_rep_3'],
            summed_col='fbf1_reads')
        p.set_sum(
            to_sum=['fbf2_rep_1', 'fbf2_rep_2', 'fbf2_rep_3'],
            summed_col='fbf2_reads')
        if re.search('fbf1', top):
            p.set_ratio(trans[dname], col2='fbf2_n2')
        elif re.search('fbf2', top):
            p.set_ratio(trans[dname], col2= 'fbf2_n2')
        else:
            print "could not tell whether %s was FBF1 or FBF2." % top
        rep_name = os.path.basename(os.path.dirname(fname))
        _mk('biological_reps_unfiltered/')
        _mk('biological_reps_filtered/')
        p.write_table('biological_reps_unfiltered/' + rep_name + '.txt')
        fp = p.get_filtered_obj(col='ratio', cutoff=cutoff)
        fp.write_table('biological_reps_filtered/' + rep_name + '.txt')

if __name__ == '__main__':
    run('fbf1_to_fbf2_n2/', cutoff=2.5)
    run('fbf2/', cutoff=2.5)
    cmd = 'python analysis/src/replicate_overlap_venn.py ' + \
' -r1 biological_reps_unfiltered/' + \
' -r2 biological_reps_unfiltered/' + \
' -o figs/biological_reps_unfiltered/'
    scribe(cmd)


