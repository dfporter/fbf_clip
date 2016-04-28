#!/usr/bin/python
import collections
import sys
import os
import re

header = ''
reads = collections.defaultdict(list)
mapqs = collections.defaultdict(int)
as_values = collections.defaultdict(int)
counted_read = set()
for li in open(sys.argv[1]).readlines():
    if li[0] == '@':
        header += li
        continue
    s = li.rstrip('\n').split('\t')
    mapq = int(s[4])
    as_flag = re.match('.*AS:i:(\d+).*', str(s[11:]))
    if as_flag is not None:
        as_value = int(as_flag.group(1))
    else:
        print "No AS value?"
    #if as_value < 10: continue
    reads[s[0]].append((mapq, as_value, s))
    if s[0] not in counted_read:
        mapqs[mapq] += 1
        as_values[as_value] += 1
        counted_read.add(s[0])
#export HISTTIMEFORMAT='%F %T  '
writ = 0
print "=" * 14
for name in reads:
    if writ > 10:
        print "Printed the first 5 reads with >1 mappping."
        break
    #if len(reads[name]) > 1:
    for i in reads[name]:
        if i[0] == 0:
            writ += 1
            print reads[name]

print "---" * 14
print 'MAPQ'
for k in mapqs:
    print "{k}:\t{o}".format(k=k, o=mapqs[k])
print "---" * 14
print 'AS'
for k in as_values:
    print "{k}:\t{o}".format(k=k, o=as_values[k])
print "Reads: %i" % len(reads)
from arrinfo import info
mapqs_l = []
for k in mapqs:
    mapqs_l.extend([float(k)] * mapqs[k])
print "~" * 14
print "MAPQ"
info(mapqs_l)
as_values_l = []
for k in as_values:
    as_values_l.extend([float(k)] * as_values[k])
print "AS: values"
info(as_values_l)
above_19 = sum([as_values[x] for x in as_values if x > 19])
below_20 = sum([as_values[x] for x in as_values if x < 20])
print "# of reads >= 20 AS:{A:>20} {p:<0}%\n# of reads < 20 AS: {B:>20}".format(
    A=above_19, B=below_20,
    p=float(100 * float(above_19)/float(above_19 + below_20)))


