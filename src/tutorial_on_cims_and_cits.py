"""

CIMS and CITS, a tutorial.

Splitting by barcode.
The first step is to split by barcode, assuming there are separate
samples in the same run, distinguished by the presence of a barcode.
In our case, we only use the forward read.
The first three bases are a random-mer.
The next four distinguish samples by barcode.

To split by barcode, we run a custom script:
    python ../clip/src/split_by_barcode.py <fastq file> <run_id> <strain>
This script will split the fastq file into four files, named:
    run_id_strain_barcode.fastq

Removing the adapter.
The 3' linker must be removed.
We call fastx_clipper to do so.
Our 3' linker is 'AGATCGGAAGAGCGGTTCAG'.
We run:
    python ../clip/src/call_fastx_clipper_to_remove_adapter.py <input_dir>
This will run the following on every *fastq file in the input_dir:
    fastx_clipper -Q 33 -n -a linker_seq -i input -o output
The ouput file will be in the directory adapter_removed, but have the sequence
    _adapter_removed.fastq
appended to the filename.
Existing filenames WILL NOT be over-written.

Checking barcode quality and trimming the barcode.


CIMS analysis.
When CIMS is called, it outputs a file with 9 columns.
Ranking by tag number (k), puts gld-1 at the top for the FBF-1 samples,
so we'll rank by tag number.
Before the output is usable, we need to:
1. Have chrm, left, right and height columns (by renaming).
    Height will be k for CIMS and score for CITS.
    Done by rename_columns().
2. Filter by FDR, and whatever else we wish to impose (k?).
    Done by apply filter.
3. Assign to a gene, and put in a column named gene_name.
    Done by assign_cims_cits_to_gene.py.
4. Include a column of the sequence named seq.
    Done by get_sequences_for_table().
5. Write a fasta.
    Done by write_fasta().
6. Write a folder with the tables named according to the analysis/src rules.
    Done by write_peaks_table()
    
We can then run the regular analysis programs.
These programs include assigning to a gene and adding the sequence.

"""
