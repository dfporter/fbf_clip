import glob
import sys
import os
import re

def reformat_collapsed_tags(
    input_dir='novo_tags_collapse', output_dir='collapsed_reformated'):
    if not os.path.exists(output_dir):
        os.system('mkdir ' + output_dir)
    for input_filename in glob.glob(input_dir + '/*'):
        print input_filename
        output_filename = output_dir + '/' + os.path.basename(input_filename)
        out_f = open(output_filename, 'w')
        with open(input_filename, 'r') as f:
            for line in f:
                s = line.rstrip('\n').split('\t')
                name = s[3].partition('[')[0]
                m = re.search('\[n=(\d+)\]', s[3])
                out_f.write("\t".join([
                    s[0], s[1], s[2], name, m.group(1), s[5]]) + '\n')

if __name__ == '__main__':
    reformat_collapsed_tags(input_dir='novo_tags_collapse/',
                            output_dir='collapsed_reformated')
