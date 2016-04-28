"""

novoalign2bed.pl:
    Convert novoalign output to bed format.
    $novoalign2bed.pl in.novoalign output.bed
    Optional parameter:
    $novoalign2bed.pl --mismatch-file fname in out
    This will set $mismatchFile to fname.
This program does the following:    
1. Counts the number of sequences in the file using wc -l.
2. For each read:
    Discards reads mapping to multiple places.
    Shift read start position by 1 up.
    Prints to $foutMismatch, if mismatchFile is set:
("\t", $chrom, $itemStart, $itemStart+1, $name,
$itemStartOnRead, $strand, $itemStartOnRead,
$refBase, ">", $mismatchBase, $matchStart), "\n" if -f $mismatchFile;
   If there are 14 columns of output, $mapDetails is set to the last column.
   $mapDetails appears to be formatted 23A>G 43T>C.
   $itemStart is the chromosomal position of these variants.
   $matchStart is apparently just 1 in all cases.
   A second format for $mapDetails occurs for deletions: 42-A 45-T
   Insertions: 36+G.
   
   For outputing tags, the following format is used:
   join(
   "\t", $chrom, $chromStart, $chromEnd, $name, $mismatch,
   $strand), "\n" if $chromStart >= 0;
   $mismatch is the number of mismatches that occur in the read.

   It may be easiest to just write our own version for sam inputs.
   sam2bed could be used for the tags, probably.

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
# My program will just call the CIMS script.

# Output files are now in ./no_barcode_plib/

python ../clip/src/call_novoalign.py

# Convert to bed format.
perl ../clip/CIMS/novoalign2bed.pl --mismatch-file mm.bed in.novo out.bed



"""
import os
import sys
import re

def sam_to_cims_beds(sam_filename, tags_filename, mismatches_filename=None):
    """Converts .sam to the two .bed files required for CIMS.
    """
    sam_file = open(sam_filename, 'r')
    tags_file = open(tags_filename, 'w')
    mismatch_file = open(mismatches_filename, 'w')
    for line_num, line in enumerate(sam_file):
        if re.match('\A@', line) is not None:
            continue  # Skip header.
        cols = line.split('\t')
        name = cols[0]
        flag = int(cols[1])
        # Flag:
        # 0 = Mapped to forward strand.
        # 16 = Mapped to reverse strand.
        # 4 = Unmapped.
        if flag == 0: strand = "+"
        elif flag == 16: strand = "-"
        else: continue
        chrom = cols[2]
        chrom_start = int(cols[3])  # 1-based leftmost.
        mapq = cols[4]
        cigar = cols[5]
        if re.search('I', cigar) is not None:
            continue
            # Skipping reads with insertions because the
            # MD tag has no clear annotation of insertions.
        num_mismatches = int(
            re.search('NM:i:(\d+)', "\t".join(cols[11:])).group(1))
        #print "num mismatches %i" % (num_mismatches)
        template_len = int(cols[8])
        seq = cols[9]
        phred_qual = cols[10]
        out_line = "\t".join([
            chrom, str(chrom_start),
            str(chrom_start + template_len),
            name, str(num_mismatches), strand])
        tags_file.write(out_line + '\n')
        parse_md(line, cols, seq, strand, mismatch_file)
        if line_num > 10000: break
    sam_file.close()
    tags_file.close()
    mismatch_file.close()

def parse_md(line, cols, seq, strand, mismatch_file):
    m = re.search("MD:Z:([\S]*)", line)
    md = m.group(1)
    parts = re.compile(
        '([0-9]*|[A-Z]|\^[A-Z]*)').split(md)
    if parts[2] != '':
        pass
        #print md, parts

    pos_in_read = 0
    for item in parts:
        if item == '': continue
        if re.match('\d+', item) is not None:
            pos_in_read += int(item)
        elif re.match('[A-Z]+', item) is not None:
            ref_seq = item
            read_seq = seq[pos_in_read:pos_in_read+len(ref_seq)]
            #print ref_seq + ">" + read_seq
            out_li = "\t".join([
                cols[2], str(pos_in_read + int(cols[3])),
                str(pos_in_read + int(cols[3]) + 1),
                cols[0], str(pos_in_read), strand,
                str(pos_in_read), ref_seq, '>', read_seq,
                '1'])
            #print out_li
            mismatch_file.write(out_li + '\n')
            pos_in_read += len(item)
        elif re.match('\^[A-Z]+', item) is not None:
            deletion = item.lstrip('\^')
            #print "Deletion of %s at pos %i" % (
            #    deletion, pos_in_read)
            out_li = "\t".join([
                cols[2], str(pos_in_read + int(cols[3])),
                str(pos_in_read + int(cols[3]) + 1),
                cols[0], str(pos_in_read), strand,
                str(pos_in_read), deletion, '-', '.',
                '1'])
            #print out_li
            mismatch_file.write(out_li + '\n')
            pos_in_read += len(deletion)
#("\t", $chrom, $itemStart, $itemStart+1, $name,
# $itemStartOnRead, $strand, $itemStartOnRead,
# $refBase, ">", $mismatchBase, $matchStart), "\n" if -f $mismatchFile;

# From plib/Bed.pm, expect bed lines:
#my @colNames = qw (
# chrom chromStart chromEnd name
# score strand thickStart thickEnd
# itemRgb blockCount blockSizes blockStarts);
# This is not the same format.

if __name__ == '__main__':
    sam_filename = 'fbf1/sams/combined_fbf1.sam'
    tags_filename = 'temp_tags.bed'
    mismatches_filename = 'temp_mismatches.bed'
    sam_to_cims_beds(
        sam_filename, tags_filename, mismatches_filename)
#    abs_path = os.path.abspath(__file__)
    path_to_plib = '/groups/Kimble/Aman\ Prasad/clip/plib/'
    os.system('export PERL5LIB=' + path_to_plib)
    os.system('perl CIMS/CIMS.pl {tag} {mismatch} {out}'.format(
        tag=tags_filename, mismatch=mismatches_filename,
        out='cims.out'))
