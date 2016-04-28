import sys
import re
import os

in_file = sys.argv[1]
outf = open(sys.argv[2], 'w')
with open(in_file, 'r') as f:
    for line in f:
        s = line.rstrip('\n').split('\t')
        if s[5] == '-':
            outf.write(line)
outf.close()
