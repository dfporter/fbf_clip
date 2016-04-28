import os
import sys
import re

seq_len_min = int(sys.argv[2])
seq_len_max = int(sys.argv[3])
out_line = ""
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if re.match('^>.*', line):
            out_line += line
        else:
            if len(line.rstrip('\n')) < seq_len_min or (
                len(line.rstrip('\n')) > seq_len_max):
                print "Error: line num %i, seq len = %i" % (
                    index,
                    len(line.rstrip('\n')))
                print "Offending line: %s" % line
                sys.exit()
            else:
                out_line += line
print out_line
