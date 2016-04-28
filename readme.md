Organization of FBF iCLIP by the original read pre-processing.
========
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57109

This repository holds the scripts and book-keeping for producing
 the FBF iCLIP paper for submission to the journal RNA.

We start with the bam files, filter, generate bed files,
 and finally generate peak calls.

We started with clip/fbf\*/bams being filtered by view -q 20 to
 clip/\*/new_filt_bams/.

We then generated the bed and bedgraph files from this, using:

```bash
sam2bed
bed_to_wig.py
```

We took RNA-seq data from clip/lib/rna_seq/*bed.

We then called the peaks using cwt-callpeaks.
We generated an auto.ini file as described in the cwt-callpeaks page.

```bash
python main.py -c auto.ini -m -n
python _load_and_combine_replicates.py  # main.py can do this.
```

Adding coverage columns, filtering and re-organizing
-------
    
```bash
cp -r fbf1/cwt-peakcaller/with_nums/* fbf1/peaks/
cp -r fbf2/cwt-peakcaller/with_nums/* fbf2/peaks/
cp -r fbf1_to_fbf2_n2/cwt-peakcaller/with_nums/* fbf1_to_fbf2_n2/peaks/
```

We want to add info columns to everything.

We're going to try making things really simple - all reads/heights will
 be normalized to per million.
 All of those values will can added by add_reads.py or
 by combine_peak_information (the latter being what was done).
 Anything not normalized will have a column heading that begins with
 unnorm_.

In the clip2/ top directory:
    
```bash
python src/combine_peak_information.py -w # Write config file
python src/combine_peak_information.py -c tmp_example_config.ini -r 3
```

This will:

1. Add maximum coverage columns in peak region in reads/million.
2. Output the filtered/ and unfiltered/ combined_fbf\*.txt files.

When run with the -w option, src/combine_peak_information.py can
 output two different config.ini files:

```python
def write_config(filename):
    # config.ini format:
    separate ='''[library]
top: /groups/Kimble/Aman Prasad/clip2/
bedgraphs_folder: %(top)s/bedgraph_unnorm/
read_beds: %(top)s/bed_collapsed/
fbf1_peak_reps: %(top)s/fbf1/peaks/
fbf2_peak_reps: %(top)s/fbf2/peaks/
fbf1_peaks: %(top)s/fbf1/peaks/combined_fbf1/null_hyp_4.txt
fbf2_peaks: %(top)s/fbf2/peaks/combined_fbf2/null_hyp_4.txt
# Bedgraph files.
fbf1_rep_1: %(top)s/fbf1/bedgraph_norm/exp_fbf1_TGGC
fbf1_rep_2: %(top)s/fbf1/bedgraph_norm/exp_fbf1_GGTT
fbf1_rep_3: %(top)s/fbf1/bedgraph_norm/exp_fbf1_CGGA
fbf2_rep_1: %(top)s/fbf2/bedgraph_norm/exp_fbf2_TGGC
fbf2_rep_2: %(top)s/fbf2/bedgraph_norm/exp_fbf2_GGTT
fbf2_rep_3: %(top)s/fbf2/bedgraph_norm/exp_fbf2_CGGA
fbf1_n2: %(top)s/fbf1/bedgraph_norm/control_n2
fbf2_n2: %(top)s/fbf2/bedgraph_norm/control_n2
rna_seq: %(top)s/bedgraph_norm/rna_seq
rna_seq_modencode: %(top)s/bedgraph_norm/rna_seq_mod4594
rna_seq_oo: %(top)s/bedgraph_norm/rna_seq_oo_srr1263137
# Bed files.
fbf1_rep_1_bed: %(top)s/fbf1/bed_collapsed/exp_fbf1_TGGC
fbf1_rep_2_bed: %(top)s/fbf1/bed_collapsed/exp_fbf1_GGTT
fbf1_rep_3_bed: %(top)s/fbf1/bed_collapsed/exp_fbf1_CGGA
fbf2_rep_1_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_TGGC
fbf2_rep_2_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_GGTT
fbf2_rep_3_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_CGGA
fbf1_n2_bed: %(top)s/bed_collapsed/fbf1_n2
fbf2_n2_bed: %(top)s/bed_collapsed/fbf2_n2
rna_seq_bed: %(top)s/bed_collapsed/rna_seq
rna_seq_modencode_bed: %(top)s/bed_collapsed/rna_seq_mod4594
rna_seq_oo_bed: %(top)s/bed_collapsed/rna_seq_oo_srr1263137
'''
    same_n2 ='''[library]
top: /groups/Kimble/Aman Prasad/clip2/
bedgraphs_folder: %(top)s/bedgraph_unnorm/
read_beds: %(top)s/bed_collapsed/
fbf1_peak_reps: %(top)s/fbf1_to_fbf2_n2/peaks/
fbf2_peak_reps: %(top)s/fbf2/peaks/
fbf1_peaks: %(top)s/fbf1_to_fbf2_n2/peaks/combined_fbf1_to_fbf2_n2/null_hyp_4.txt
fbf2_peaks: %(top)s/fbf2/peaks/combined_fbf2/null_hyp_4.txt
# Bedgraph files.
fbf1_rep_1: %(top)s/fbf1_to_fbf2_n2/bedgraph_norm/exp_fbf1_TGGC
fbf1_rep_2: %(top)s/fbf1_to_fbf2_n2/bedgraph_norm/exp_fbf1_GGTT
fbf1_rep_3: %(top)s/fbf1_to_fbf2_n2/bedgraph_norm/exp_fbf1_CGGA
fbf2_rep_1: %(top)s/fbf2/bedgraph_norm/exp_fbf2_TGGC
fbf2_rep_2: %(top)s/fbf2/bedgraph_norm/exp_fbf2_GGTT
fbf2_rep_3: %(top)s/fbf2/bedgraph_norm/exp_fbf2_CGGA
fbf1_n2: %(top)s/fbf1_to_fbf2_n2/bedgraph_norm/control_n2
fbf2_n2: %(top)s/fbf2/bedgraph_norm/control_n2
rna_seq: %(top)s/bedgraph_norm/rna_seq
rna_seq_modencode: %(top)s/bedgraph_norm/rna_seq_mod4594
rna_seq_oo: %(top)s/bedgraph_norm/rna_seq_oo_srr1263137
# Bed files.
fbf1_rep_1_bed: %(top)s/fbf1_to_fbf2_n2/bed_collapsed/exp_fbf1_TGGC
fbf1_rep_2_bed: %(top)s/fbf1_to_fbf2_n2/bed_collapsed/exp_fbf1_GGTT
fbf1_rep_3_bed: %(top)s/fbf1_to_fbf2_n2/bed_collapsed/exp_fbf1_CGGA
fbf2_rep_1_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_TGGC
fbf2_rep_2_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_GGTT
fbf2_rep_3_bed: %(top)s/fbf2/bed_collapsed/exp_fbf2_CGGA
fbf1_n2_bed: %(top)s/bed_collapsed/fbf1_n2
fbf2_n2_bed: %(top)s/bed_collapsed/fbf2_n2
rna_seq_bed: %(top)s/bed_collapsed/rna_seq
rna_seq_modencode_bed: %(top)s/bed_collapsed/rna_seq_mod4594
rna_seq_oo_bed: %(top)s/bed_collapsed/rna_seq_oo_srr1263137
'''
    open(filename, 'w').write(separate)
```

Adding annotations
-------

Then add the needed annotation information:
    
```bash
python analysis/src/apply_annotatePeaks.py -r -c tmp_example_config.ini -i DIR
# Given -r, this adds "ratio" as a ratio of raw read number,
# which is what we want for test_various_controls and the 5/6 list.
# Run again after that to get the normalized numbers used for everything
# else.
nohup python analysis/src/apply_annotatePeaks.py -r -c tmp_example_config.ini -i fbf1/peaks/ > qq 2>&1 </dev/null & 
nohup python analysis/src/apply_annotatePeaks.py -r -c tmp_example_config.ini -i fbf2/peaks/ > qq 2>&1 </dev/null & 
nohup python analysis/src/apply_annotatePeaks.py -r -c tmp_example_config.ini -i fbf1_to_fbf2_n2/peaks/ > q3 2>&1 </dev/null &
cp -r fbf1/peaks/* fbf1/backup_of_peaks/
cp -r fbf2/peaks/* fbf2/backup_of_peaks/
cp -r fbf1_to_fbf2_n2/peaks/* fbf1_to_fbf2_n2/backup_of_peaks/
```

This will:

    1. Add location, seq, has_fbe, biotype.
    2. Add sum (actually the average) and ratio columns.

Make the 5/6 list
-------
```bash
# This moves things into the for_methods_comparison/ directory.
# It also creates the 5/6 filtered files.
python analysis/src/test_various_controls.py -g config.ini -b
cp for_methods_comparison/filtered/4_five_reps/combined_fbf.txt filtered_five_reps/
# Make the figures.
python analysis/src/test_various_controls.py -g config.ini -c
```

Make the filtered and unfiltered FBF-1/-2 lists
-------

To make the separate lists, we edit combine_peak_information.py to write
 the other config.ini file, and then run:

```bash
python src/combine_peak_information.py -w # Write config file
python src/combine_peak_information.py -c tmp_example_config.ini -r 3
# No -r operator, so the sums and ratios and now for normalized reads.
python analysis/src/apply_annotatePeaks.py -c tmp_example_config.ini -i unfiltered/
python analysis/src/apply_annotatePeaks.py -c tmp_example_config.ini -i filtered/
# This small script makes sure the fbf1_n2 reads columns in *filtered/combined_fbf1.txt
# are from the normalized FBF-1 N2 sample.
python src/add_fbf1_n2_reads_column.py
# This is needed to add the RNA-seq columns for Table S2.
python src/add_rna_seq_data.py

```


Making figure 2
-------
make_figs.py now tries to annotate everything itself, and can act on the
 output of combine_peak_information.py. Edit the code before using.

Calling `test_various_controls.py -g config.ini -b` has the following effect:

1. For each null_hyp peaks file in `lib['fbf*_rep_*_all_hyp']` from
the config.ini file, copy it to for_methods_comparison/separate_reps/.
2. Call analysis/src/subset_targs_by_all_replicates.py on each of these
folders.
3. If 'seq' or 'has_fbe' columns are missing, add info to the resulting 
/for_methods_comparison/unfiltered/ folders by calling
 add_info_columns.py -e -s -d {input_dir} on each.
4. Create the filtered/ folders by using a hard-coded cutoff to the ratio
column.

Info for text
-------

```bash
# In top directory clip2/
python analysis/src/info_for_text_ui.py -f filtered_five_reps/ -s filtered/
```

Supplementary tables
-------

Should be included with make_figs.py?

```bash
python analysis/src/create_xls_tables.py
python analysis/src/make_table_s4.py
```

Bernstein figure
-------
```bash
python analysis/src/feature_locations/mutations_in_fbe.py IN_FASTA
# Writes to figs/Bernstein_*.pdf

# Used:
python analysis/src/feature_locations/mutations_in_fbe.py data/fasta/5reps_subpeaks/highest_both.fa 
```

Score metrics
-------
```bash
# In cwt directory.
python score_metrics.py -c auto.ini -p ../../filtered/

With ncRNA:
    ---
Dataset: combined_fbf1.txt
Number of peaks: 2946
Number of genes: 2345
With FBE: 2043, 69.0%
Without FBE: 903
Positive controls: 13/15
Missing positives: set(['fog-3', 'egl-4'])

Dataset: combined_fbf2.txt
Number of peaks: 1564
Number of genes: 1457
With FBE: 1319, 84.0%
Without FBE: 245
Positive controls: 13/15
Missing positives: set(['fog-3', 'egl-4'])

Without ncRNA:
    ---
Dataset: combined_fbf1.txt
Number of peaks: 2878
Number of genes: 2290
With FBE: 2002, 69.0%
Without FBE: 876
Positive controls: 13/15
Missing positives: set(['fog-3', 'egl-4'])

Dataset: combined_fbf2.txt
Number of peaks: 1530
Number of genes: 1424
With FBE: 1292, 84.0%
Without FBE: 238
Positive controls: 13/15
Missing positives: set(['fog-3', 'egl-4'])

with ncRNA:
---
Dataset: combined_fbf.txt
Number of peaks: 1609
Number of genes: 1404
With FBE: 1252, 77.0%
Without FBE: 357
Positive controls: 13/15
Missing positives: set(['fog-3', 'egl-4'])

without ncRNA:
---
python score_metrics.py -c auto.ini -p ../../filtered_five_reps/
Dataset: combined_fbf.txt
Number of peaks: 1569
Number of genes: 1369
With FBE: 1228, 78.0%
Without FBE: 341
Positive controls: 13/15
Missing positives: set(['fog-3', 'egl-4'])


# make_figs.py adds biotype, location, ect.
python make_figs.py filtered/
python make_figs.py unfiltered/

```

Alternative to combine_peak_inforamation:
```
# In cwt directory.
python ../../src/add_reads.py -c auto.ini
python filter.py -p fbf2/peaks/combined_fbf2/null_hyp_4.txt -c auto.ini -r 10 -o filtered/fbf2.txt

# Add fbf\*_reads and n2_from_fbf\*_reads.
```

Ortholist comparison
-------

Our strategy is to subset the worm and human genomes to those
worm loci with Ortholist ENSGxx and human ENSGxx with Ortholist
worm loci. So a complexTranslation object is created holding
just those loci and ENSGxx. We then pare down the FBF and
PUM2 target lists (as complexTargetSets) to those loci and
ENSGxx that exist in Ortholist. We can then get the overlap
information with overlap_with_complexTargetSet.

This code is all contained in ortho/.   

```bash
python pum2_vs_fbf.py
python yeast_vs_fbf.py
```




