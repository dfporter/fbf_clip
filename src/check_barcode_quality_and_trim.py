"""
CIMS processing pathway:
Run in the folder containing the raw fastq files in subfolders:
python ../clip/src/split_by_barcode.py

python ../clip/src/call_fastx_clipper_to_remove_adapter.py

# This reuquires the first 9 bases each have 20 quality, then removes them.
python ../clip/src/check_barcode_quality_and_trim.py
# The CIMS equivalent program is stripBarcode.pl,
# which can be run as:
export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/
perl ../clip/CIMS/stripBarcode.pl -len 9 -format fastq in.fastq out.fastq
# This program will filter by the barcode quality scores
# and then just call the CIMS script.

# Output files are now in ./no_barcode_plib/

python ../clip/src/call_novoalign.py

# Convert to bed format.
perl ../clip/CIMS/novoalign2bed.pl --mismatch-file mm.bed in.novo out.bed


"""
import glob
import sys
import os
import re
import HTSeq

# Require minimum of 20 base quality at barcode and
# random-mer positions, or the first 3 + 4 + 2 = 9 bases.
# Previously, used fastx_trimmed -f 10.
# We'll check the quality of these bases first.

# Also want to filter by total MAPQ 20 later, if using bowtie.


os.system('mkdir no_adapter_no_barcode')
os.system('mkdir no_barcode_plib')
cmd = r'''export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/'''
print cmd
os.system(cmd)
os.system('mkdir barcode_qual_filter/')
for input_fastq in glob.glob('adapter_removed/run813*.fastq'):
    print input_fastq
    output_filename = 'no_barcode_plib/' + os.path.basename(input_fastq).partition('.fastq')[0]
    output_filename += '_no_barcode.fastq'
    if os.path.exists(output_filename):
        continue
    use_plib = True
    filter_barcode_qual = True
    if use_plib:
        if filter_barcode_qual:
            filtered_filename = 'filtered_novo_tags/'
            filtered_filename += os.path.basename(input_fastq)
            if not os.path.exists(filtered_filename):
                fastq_out = open(filtered_filename, 'w')
                fastq = HTSeq.FastqReader(input_fastq)
                for read in fastq:
                    if min(read.qual[0:10]) < 20: continue
                    else: read.write_to_fastq_file(fastq_out)
            if os.path.exists(output_filename): continue
            cmd = 'perl ../clip/CIMS/stripBarcode.pl -len 9 -format fastq'
            cmd += ' %s %s' % (filtered_filename, output_filename)
            print cmd
            os.system(cmd)            
        else:    
            cmd = 'perl ../clip/CIMS/stripBarcode.pl -len 9 -format fastq'
            cmd += ' %s %s' % (input_fastq, output_filename)
            print cmd
            os.system(cmd)
    else:
        fastq_out = open(output_filename, 'w')
        fastq = HTSeq.FastqReader(input_fastq)
        for read in fastq:
            # read is a SequenceWithQualities object.
            if min(read.qual[0:10]) < 20:
                continue
            else:
                # Trim and write out.
                read_out = HTSeq.SequenceWithQualities(
                    read.seq[10:],
                    read.name,
                    read.qualstr[10:])
    #            read.seq = read.seq[10:]
    #            read.qualstr = read.qualstr[10:]
    #            read.qual = read.qual[10:]
                read_out.write_to_fastq_file(fastq_out)
        fastq_out.close()
