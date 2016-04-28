import pandas
import numpy as np

import scipy as sp
import sys
import glob
indir = sys.argv[1]

print ''
print ''
print 'Stats:'
for fname in glob.glob('methods/*/comb*'):
    print '*' * 7
    print fname
    #try:
    df = pandas.read_csv(fname, sep='\t')
    stats = df['neg_ip_local_norm'].tolist()
    stats = np.array(stats)
    #fbf1_reads = df['fbf1_reads'].tolist()
    #fbf2_reads = df['fbf2_reads'].tolist()
    print "Minimum p value {m}. Maximum p value {ma}.".format(
        m=min(stats), ma=max(stats)
    )
    #print "Minimum fbf1_reads {mf}. Minimum fbf2_reads {mf2}.".format(
    #    mf=min(fbf1_reads), mf2=min(fbf2_reads)
    #)
    #except:
    #    print 'failed on ' + fname
