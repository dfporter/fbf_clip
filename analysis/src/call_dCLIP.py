import os
import glob
import re
import sys

top_level_dir = sys.argv[1]
out_dir = 'data/' + os.path.basename(top_level_dir.rstrip('/'))
if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
out_dir += '/dCLIP_out'
if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
cmd = 'perl dCLIP1.7/bin/dCLIP.pl -f1 {} -f2 {} -m1 100 -m2 100'.format(
    'fbf1/sams/combined_fbf1.sam', 'fbf2/sams/combined_fbf2.sam')
cmd += ' -dir {} -iCLIP 2'.format(out_dir)
print cmd
os.system(cmd)
