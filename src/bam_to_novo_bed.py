__author__ = 'dp'
import HTSeq
import os
import sys
import re
import glob
import argparse


def move_barcode_to_name_in_fastq(filename, out_dir):
    if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
    outf = open(out_dir + '/%s' % os.path.basename(filename), 'w')
    fastq = HTSeq.FastqReader(filename)
    obs_let = set()
    # phred 33
    for read in fastq:
        if len(read.seq) < 14: continue
        if min(read.qualstr[:9]) < 53: continue
        n_read = HTSeq.SequenceWithQualities(
            read.seq[9:],
            read.name.partition(' ')[0] + '#' + read.seq[0:9],
            read.qualstr[9:])
        n_read.write_to_fastq_file(outf)
    

def convert_sam_to_novo_bed(filename, out_dir):
    if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
    outf = open(
        out_dir + '/%s' % os.path.basename(filename).partition('.sam')[0] + '.bed', 'w')
    with open(filename, 'r') as f:
        for read in f:
            if read[0] == '@': continue
            s = read.rstrip('\n').split('\t')
            if s[1] == '4': continue
            chrom = s[2]
            left = int(s[3])
            #mapq = s[4]
            if s[1] == '16':
                strand = '-'
                #print read
            else:
                strand = '+'
            right = str(left + len(s[9]))
            if int(left) > int(right):
                print read
            num_mut = 0
            cigar = s[5]
            insertions = 0
            for mut in re.finditer('(\d+)I', s[5]):
                print mut.groups()[0]
                insertions += int(mut.groups()[0])
            right -= insertions
            for mut in re.finditer('[ID]', s[5]): num_mut += 1
            outli = "{chrm}\t{l}\t{r}\t{name}\t{score}\t{strand}\n".format(
                chrm=s[2], l=s[3], r=right, name=s[0], score=num_mut, strand=strand)
            outf.write(outli)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Move barcodes or convert sam to bed for CIMS.''')
    parser.add_argument('-m', '--move_barcode',
                     default=False, action='store_true',
                     help='Move barcodes in fastq files in the input directory to the read names.',
                    )
    parser.add_argument('-s', '--sam',
                     default=False, action='store_true',
                     help='Convert sam to .bed for CIMS')
    parser.add_argument('-i', '--input_dir',
                     help='Input directory.')
    args = parser.parse_args()
    if args.move_barcode:
        for filename in glob.glob(args.input_dir + '/*.fastq'):
            move_barcode_to_name_in_fastq(filename, 'adapter_moved_to_name/')
    elif args.sam:
        for filename in glob.glob(args.input_dir + '/*.sam'):
            convert_sam_to_novo_bed(filename, 'bed/')
