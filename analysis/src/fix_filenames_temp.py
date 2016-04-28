import re
import sys
import glob
import os

indir = sys.argv[1]

for filename in glob.glob(indir + '/*'):
    s = filename.partition('combined')
    out_base = "".join([s[1], s[2]])
    out_name = s[0] + '/' + out_base
    print out_name
    cmd = 'cp {i} {d}/{on}'.format(
        i=filename, d=os.path.split(s[0])[1],
        on=out_base)
    print cmd
    os.system(cmd)
