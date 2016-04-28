import glob
import os
import sys


#input_dir = sys.argv[1]

def convert_folder(input_dir):
    for filename in glob.glob(input_dir + '/*.bam'):
        print filename
        out_name = os.path.dirname(os.path.dirname(input_dir)) + '/'
        out_name += 'beds/'
        out_name += os.path.basename(filename).partition('.bam')[0] + '.bed'
        cmd = 'bamToBed -i %s > %s' % (filename, out_name)
        print cmd
        os.system(cmd)

convert_folder('fbf1/bams/')
convert_folder('fbf2/bams/')
