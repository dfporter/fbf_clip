import glob
import sys
import os
import re

l3 = 'AGATCGGAAGAGCGGTTCAG'
input_dir = sys.argv[1]
input_filenames = glob.glob('%s/*fastq' % input_dir)
if not os.path.exists('adapter_removed/'):
    os.system('mkdir adapter_removed')
for input_filename in input_filenames:
    output_filename = 'adapter_removed/' + os.path.basename(
        input_filename).partition('.fastq')[0] + '_adapter_removed.fastq'
    if os.path.exists(output_filename):
        continue
    if re.search('_adapter_removed.fastq', input_filename) is not None:
        continue
    cmd = 'fastx_clipper -Q 33 -n -a %s -i %s -o %s' % (
        l3, input_filename, output_filename)
    print cmd
    os.system(cmd)    
