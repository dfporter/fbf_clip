import glob
import sys
import os
import re

if not os.path.exists('novo_tags'):
    os.system('mkdir novo_tags')
if not os.path.exists('mismatch_tags'):
    os.system('mkdir mismatch_tags')
for input_filename in glob.glob('novoaligned/*.novo'):
    tags_filename = 'novo_tags/' + os.path.basename(
        input_filename).partition('.novo')[0] + '.bed'
    mismatch_filename = 'mismatch_tags/' + os.path.basename(
        input_filename).partition('.novo')[0] + '.bed'
    cmd = '../clip/CIMS/novoalign2bed.pl --mismatch-file {mm} {i} {o}'.format(
        mm=mismatch_filename, i=input_filename, o=tags_filename)
    print cmd
    os.system(cmd)
