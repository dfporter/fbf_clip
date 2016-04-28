import glob
import sys
import os
import re

os.system('mkdir novo_tags_collapse')
#for input_filename in glob.glob('novo_tags/*'):
for input_filename in glob.glob('bed/*'):
    output_filename = 'novo_tags_collapse/'
    output_filename += os.path.basename(input_filename)
    if os.path.exists(output_filename): continue
    cmd = 'perl ../../Aman\ Prasad/clip/CIMS/tag2collapse.pl --random-barcode --seq-error-model fix=0.02 -EM 50 --keep-cache {i} {o}'.format(
        i=input_filename, o=output_filename)
    #cmd = 'perl ../clip/CIMS/tag2collapse.pl --random-barcode --seq-error-model fix=0.02 -EM 50 --keep-cache {i} {o}'.format(
    #    i=input_filename, o=output_filename)
    print cmd
    os.system(cmd)
