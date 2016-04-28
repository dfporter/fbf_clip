"""

Run in the folder containing the raw fastq files in subfolders:
python ../clip/src/split_by_barcode.py

python ../clip/src/call_fastx_clipper_to_remove_adapter.py

# This reuquires the first 9 bases each have 20 quality, then removes them.
python ../clip/src/check_barcode_quality_and_trim.py

# Output files are now in ./no_barcode_plib/

python ../clip/src/call_novoalign.py

#Output files are now in ./novoaligned/.

perl ../clip/CIMS/novoalign2bed.pl --mismatch-file mm.bed novoaligned/in.novo out.bed

"""

import glob
import sys
import os
import re

input_filenames = glob.glob('no_barcode_plib/*.fastq')

for input_filename in input_filenames:
    output_filename = './novoaligned/' + os.path.basename(
        input_filename).partition('.fastq')[0] + '.novo'
    if os.path.exists(output_filename):
        print "call_novoalign(): Output file %s already exists. Skipping..." % (
            output_filename)
        continue
    cmd = '/opt/novocraft/novoalign -d /scratch/indexes/ws235.novo'
    cmd += ' -f {infile} > {outfile}'.format(
        infile=input_filename, outfile=output_filename)
    print cmd
    os.system(cmd)
