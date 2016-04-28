#!/usr/bin/env bash
# Input argument: name of the top level directory where all
# the peaks are. Example: ./pypeaks_fdr5_negip_local

# This document is intended to function, in its final version,
# as a bash script that will take any set of called peaks,
# (and stable library files), and run everything required to
# do all analytics and produce all figures.
# The reason for this is the high probabilty that peak lists
# will change over the course of paper writing, and it will be
# impossible to make progress if everything needs to be re-run
# by hand when something at the start of the pipeline changes.

# We'll place analytical src code in analytics/src/.

# Peak files need to be ordered such that a different top level
# directory can be specified as a single input to each analytical
# program, and all of the correct peaks files for that method
# of peak calling are automatically found.
# We'll do this by taking putting each method of calling peaks
# under a peaks/method/ folder. FBF-1 and FBF-2 are then separated into
# peaks/method/combined_fbf1.txt and peaks/method/combined_fbf2.txt files, for
# peaks from combined replicates. Individual replicates are
# put in the same folder, but are given the specific names
# peaks/method/fbf1_aacc.txt, or the analogous name.

# To take the peaks from fbf*/peaks/ folders into a single methods folder,
# run python callpeaks_v2/create_dirs_for_analysis.py.

method=$1

case "$method" in
*/)
    wslash=$1
    ;;
*)
    wslash=$1/
    ;;
esac


method=$1




# Not in a figure:
#     A. Add info columns.
#     B. Create a list of 5/6 replicates.
#     C. Compare whether multiple peaks occur more in longer transcripts.
#     D. Add bed info.
#     E. Split by reads-in-peak filter.

# A. Add info columns.
echo
#echo python analysis/src/add_info_columns.py -d ${method} -e -w -c
#python analysis/src/add_info_columns.py -d ${method} -e -w -c

# B. 5/6 replicates.
# I note add_info.py also does this if given the -c parameter. Are they the same?
echo
#echo python analysis/src/subset_targs_by_all_replicates.py ${method}
#python analysis/src/subset_targs_by_all_replicates.py ${method}

# Add info columns to the 5/6 replicates separate column.
echo
#echo python analysis/src/add_info_columns.py -d ${method%/}_five_reps/ -e -w -c
#python analysis/src/add_info_columns.py -d ${method%/}_five_reps/ -e -w -c

# C. Do peak pairs occur more in longer transcripts?
filedir=${wslash%/*}
if
echo $filedir | grep five
then
var=${filedir}/combined_fbf.txt
else
var=${filedir%/}_five_reps/combined_fbf.txt
fi
echo
#echo python analysis/src/peak_pair_locations_ui.py --no_user_input -p $var
#python analysis/src/peak_pair_locations_ui.py --no_user_input -p $var


# E. Add reads-in-peak numbers.
# filter and splitting the peaks that pass the filter to a new folder.
# Use -d to keep from over-writting columns.
echo
#echo python src/filter_peaks_by_ratio_in_peak_region.py -d --load_bed_info -i ${method} --output_dir ${method}
#python src/filter_peaks_by_ratio_in_peak_region.py -d --load_bed_info -i ${method} --output_dir ${method}

case "$method" in
*filtered*/)
    is_filtered=1
    ;;
*)
    is_filtered=0
    # E. Split the input folder by applying a reads-in-peak ratio
    echo
    echo python src/filter_peaks_by_ratio_in_peak_region.py --analysis -i ${method} --output_dir filtered_${method}  --use fbf2_n2
    python src/filter_peaks_by_ratio_in_peak_region.py --analysis -i ${method} --output_dir filtered_${method} --use fbf2_n2
    ;;
esac


# Fig 2: Basic iCLIP data quality metrics:
#    A. Table of % unique reads mapped by replicate.
#    B. Replicate overlap.
#    C. Scatterplot comparison to RNA-seq.
#    D. UGUNNNAU and -1/-2C FBEs are enriched in the highest peaks.
#    E. Venn of overlap with RIP-ChIP.

# Fig 2.A: Table of unique reads and mapping statistics.
echo
#echo python analysis/src/table_of_dataset_sizes_and_statistics.py
#python analysis/src/table_of_dataset_sizes_and_statistics.py
# Outputs to tables/mapping_quality_stats.xls.

# Fig 2.B gene names for Fig 3.C:
echo
echo python analysis/src/replicate_overlap_venn.py ${method}
python analysis/src/replicate_overlap_venn.py ${method}
# replicate_overlap_venn.py also produces the FBF-1/2 overlap venn.
# It also produces the table of unique and separate targets.

# Fig 2.C: Scatterplot vs RNA-seq.
echo
#echo python analysis/src/scatterplot_correlation_with_controls.py --cols --input ${method}
#python analysis/src/scatterplot_correlation_with_controls.py --cols --input ${method}

# Fig 2.D: FBE vs RANK.
echo
echo python analysis/src/plot_fbe_vs_rank.py ${method}
python analysis/src/plot_fbe_vs_rank.py ${method}
# Outpus to figs/method/basename_one_to_two.eps.

# Fig 2.E Overlap with RIP-ChIP.
echo
echo python analysis/src/venn_of_comparison_to_rip_chip.py ${method}
python analysis/src/venn_of_comparison_to_rip_chip.py ${method}
# Writes to figs/${method}/overlap_with_rip_*.

# Fig 3: FBF-1 vs FBF-2.
#     A. FBF-1/FBF-2 scatterplot correlation.
#     B. FBF-1/FBF-2 correlation heatmap.
#     C. patr-1 or other RNA showing differences.

# Fig 3.A: Scatterplot of correlation between FBF-1 and FBF-2.
# -p option here will load existing data to save time.
echo
#echo python analysis/src/scatterplot_correlation_by_wig_ui.py --no_ui -d ${method}
#python analysis/src/scatterplot_correlation_by_wig_ui.py --no_ui -d ${method}
echo
#echo python analysis/src/scatterplot_correlation_fbf1_and_fbf2.py -d ${method}
#python analysis/src/scatterplot_correlation_fbf1_and_fbf2.py -d ${method}
# Outputs to /figs/method/scatterplot_by_wig.pdf and
# data/method/outliers_fbf1_vs_fbf2.txt.

# Fig 3.B: Tables are output with Fig 2.2, and need GO run on them.
echo
#echo python analysis/src/extract_ranges_and_write_as_bed.py ${method}
#echo skipping...
#python analysis/src/extract_ranges_and_write_as_bed.py ${method}
echo
#echo python analysis/src/call_bamCorrelate.py ${method}
#echo skipping...
#python analysis/src/extract_ranges_and_write_as_bed.py ${method}
echo
#python analysis/src/call_bamCorrelate.py ${method}
# Outputs to figs/method/spearman_correlation_between_bams_in_peak_regions.png.

# Fig 3.C: Get from IGV.

# Fig 3.4: Tables of RNA types.
echo
echo python analysis/src/table_of_rna_types.py ${method}
python analysis/src/table_of_rna_types.py ${method}

# Fig 4: Pooled FBF data reveal high-confidence targets.
#    A. Pathway analysis.
#    B. Peak traces.
#    C. Subset of conserved targets.

# Fig 5: How FBF-1 and FBF-2 bind their targets.
#    A. CIMS/CITS cartoon of method.
#    B. XL sites relative to FBE.
#    C. Motif enrichments by MEME.
#    D. Motif enrichments from the peak caller.
#    E. Alternate FBEs.

# Fig 5.A: Made by hand.

# Fig 5.B: XL sites relative to FBE.
# Have uncomment everything needed from the masive cims_master.py
# pipeline. cims_master.py will also call
# cims_pos_vs_fbe.cims_vs_site_in_dir on the folder
# reproducible/. Running this stuff takes some attention paid
# to what lines should be uncommented.
echo Running CIMS and CITS analysis...
echo skpping
#export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/
#(cd ../raw\ fastq\ for\ Doug\ CIMS\ CITS/ && exec python ../src/cims_master.py)
#(cd ../raw\ fastq\ for\ Doug\ CIMS\ CITS/ && exec python ../src/cims_pos_vs_fbe.py)

# Fig 5.C: MEME-ChIP results.
# Run on the reproducible CIMS/CITS .fasta files put in
# raw\ fastq\ for\ Doug\ CIMS\ CITS/fasta/.

# Fig 5.D: MEME on the peak caller output.
# Fasta files are in fasta/.

# Fig 5.E: Alternative FBEs.
echo
#echo python analysis/src/alternative_fbes.py ${method}
#echo skipping...
#python analysis/src/alternative_fbes.py ${method}
# Outputs to data/method/cryptic_fbes_numofshuffles.txt.

# Fig 6: Binding site location.
#     A. Pie chart of binding location.
#     B. Position in a normalized transcript.
#     C. Peak traces for 5 candidate 5 prime UTR targets.
#     D. Peak location in 3 prime UTRs.

# Fig 6.A: Pie chart of binding site location.
echo
echo python analysis/src/peak_locations.py ${method}
python analysis/src/peak_locations.py ${method}
# Outputs to figs/method/peak_locations_pie_chart_label.pdf.

# Fig 6.B, 6.D: Position in a normalized transcript.
echo
echo python analysis/src/feature_locations/determine_feature_locations.py -i $var
python analysis/src/feature_locations/determine_feature_locations.py -i $var
# Outputs to figs/features_in_normalized_gene.pdf,
# figs/features_in_utr.pdf


# Fig 6.C: Get from IGV.

# Fig 7: Different methods of controlling iCLIP data.
#     A. Three barplots for the methods of initial control.
#     B. Barplot of one method vs a gradient of enrichment ratio.

# Fig 7.A: Different methods of initial control.
echo
#echo python analysis/src/test_various_controls.py --build_from_raw_peaks
#python analysis/src/test_various_controls.py --build_from_raw_peaks
#echo python analysis/src/test_various_controls.py --info
#python analysis/src/test_various_controls.py --info
#echo python analysis/src/test_various_controls.py --combined_fig
#python analysis/src/test_various_controls.py --combined_fig
# Outputs to figs/various_controls/effect_of_controls_barplot_.pdf

# Fig 7.B: Gradient of entichment vs cutoff.
echo
#echo python analysis/src/effect_of_different_controls_figure.py --input ${method}/combined_fbf.txt
#echo skipping...
#python analysis/src/effect_of_different_controls_figure.py --input ${method}/combined_fbf.txt
# Outputs to figs/method/effect_of_controls_barplot_*.pdf


# Currently unused: dCLIP.
# dCLIP needs SAM files.
echo
#echo python analysis/src/call_dCLIP.py ${method}
#python analysis/src/call_dCLIP.py ${method}
#echo python analysis/src/get_fasta_of_dCLIP-run_DREME.py ${method}
#python analysis/src/get_fasta_of_dCLIP-run_DREME.py ${method}



