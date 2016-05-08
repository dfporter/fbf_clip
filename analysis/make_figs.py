
import sys
import os
import re
import glob
import pandas
import HTSeq

sys.path.insert(
    0, '/groups/Kimble/Aman Prasad/clip2/src/')
sys.path.insert(
    0, '/groups/Kimble/Aman Prasad/clip2/analysis/src/')
import biological_replicates_folder 
import apply_annotatePeaks
del sys.path[0]

def make_figs(method, has_info=True):
    # Get data for annotation.
    #gtf = pandas.read_csv('lib/gtf_with_names_column.txt', sep='\t')
    #fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    #sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    #gtf_r = gtf.to_dict('records')
    #gtf_d = collections.defaultdict(list)
    #for row in gtf_r:
    #    gtf_d[row['gene_name']].append(row)
    #ga = load_bedgraph('bedgraph_norm/combined_fbf1.wig')
    # Determine what kind of folder this is.
    with_slash = method.rstrip('/') + '/'
    without_slash = method.rstrip('/')
    is_filtered = False
    is_five_reps = False
    if re.match('.*filtered.*', method) is not None:
        is_filtered = True
    if re.match('.*five_rep.*', method) is not None:
        is_five_reps = True
    #for fname in glob.glob(with_slash + '*.txt'):
    #    apply_annotatePeaks.add_info(fname, gtf, sequences, ga, gtf_d)
    if not os.path.exists('figs/%s/' % os.path.basename(os.path.dirname(with_slash))):
        os.system('figs/%s/' % os.path.basename(os.path.dirname(with_slash)))
    # Not in a figure:
    #     A. Add info columns.
    #     B. Create a list of 5/6 replicates.
    #     C. Compare whether multiple peaks occur more in longer transcripts.
    #     D. Add bed info.
    #     E. Split by reads-in-peak filter.

    for fname in glob.glob(with_slash + '*.txt'):
        peaks = pandas.read_csv(fname, sep='\t')
        if is_five_reps:
            msg = """
Add height to %s as average of fbf1_reads and fbf2_reads? 
This will be an average of the individual reads per million values.
(Y/N):""" % fname
            if raw_input(msg).upper() == 'Y':
                tups = zip(peaks['fbf1_reads'].tolist(), peaks['fbf2_reads'].tolist())
                peaks['height'] = [sum(list(x))/2.0 for x in tups]
                peaks.to_csv(fname, sep='\t', index=False)
        else:
            if re.search('fbf1', fname) and not (re.search('fbf2', fname)):
                msg = "Add height column to %s as fbf1_reads? (Y/N):" % fname
                if raw_input(msg).upper() == 'Y':
                    peaks['height'] = peaks['fbf1_reads']
                    peaks.to_csv(fname, sep='\t', index=False)
            if re.search('fbf2', fname) and not (re.search('fbf1', fname)):
                msg = "Add height column to %s as fbf2_reads? (Y/N):" % fname
                if raw_input(msg).upper() == 'Y':
                    peaks['height'] = peaks['fbf2_reads']
                    peaks.to_csv(fname, sep='\t', index=False)
    if (not has_info) and (not is_filtered):
        # A. Add info columns.
        scribe('python analysis/src/add_info_columns.py -d %s -e -w' % with_slash)
    #if not is_five_reps:
        # B. 5/6 replicates.
        # I note add_info.py also does this if given the -c parameter. Are they the same?
        #scribe('python analysis/src/subset_targs_by_all_replicates.py %s' % with_slash)

        # Add info columns to the 5/6 replicates separate column.
        #scribe('python analysis/src/add_info_columns.py -d %s_five_reps/ -e -w' % without_slash)

    # C. Do peak pairs occur more in longer transcripts?
    scribe('python analysis/src/peak_pair_locations_ui.py --no_user_input -p %s' % with_slash)
    if not is_filtered:
        # E. Add reads-in-peak numbers.
        # filter and splitting the peaks that pass the filter to a new folder.
        # Use -d to keep from over-writting columns.
        scribe('python src/filter_peaks_by_ratio_in_peak_region.py -d --load_bed_info -i %s --output_dir %s' % (
            with_slash, with_slash))
        scribe('python src/filter_peaks_by_ratio_in_peak_region.py \
--analysis -i %s --output_dir filtered_%s --use fbf2_n2' % (
            with_slash, with_slash))
    # Fig 2: Basic iCLIP data quality metrics:
    #    A. Table of % unique reads mapped by replicate.
    #    B. Replicate overlap.
    #    C. Scatterplot comparison to RNA-seq.
    #    D. UGUNNNAU and -1/-2C FBEs are enriched in the highest peaks.
    #    E. Venn of overlap with RIP-ChIP.

    # Fig 2.A: Table of unique reads and mapping statistics.
    #scribe('python analysis/src/table_of_dataset_sizes_and_statistics.py')
    # Outputs to tables/mapping_quality_stats.xls.
    if not is_five_reps:
        # Fig 2.B gene names for Fig 3.C:
        os.system('echo')
        biological_replicates_folder.run('fbf1_to_fbf2_n2/', cutoff=3)
        biological_replicates_folder.run('fbf2/', cutoff=3)
        #scribe('python analysis/src/replicate_overlap_venn.py -c %s -o figs/%s' % (
        #    with_slash, os.path.basename(os.path.dirname(with_slash))))
        if not os.path.exists('figs/biological_reps_unfiltered/'):
            os.system('figs/biological_reps_unfiltered/')
        if not os.path.exists('figs/biological_reps_filtered/'):
            os.system('figs/biological_reps_filtered/')
        cmd = 'python analysis/src/replicate_overlap_venn.py ' + \
' -r1 biological_reps_unfiltered/' + \
' -r2 biological_reps_unfiltered/' + \
' -o figs/biological_reps_unfiltered/' 
        scribe(cmd)
#        cmd = 'python analysis/src/replicate_overlap_venn.py ' + \
#' -r1 biological_reps_filtered/' + \
#' -r2 biological_reps_filtered/' + \
#' -o figs/biological_reps_filtered/' 
#        scribe(cmd)
        # replicate_overlap_venn.py also produces the FBF-1/2 overlap venn.
        # It also produces the table of unique and separate targets.
        # Fig 2.C: Scatterplot vs RNA-seq.
        os.system('echo')
        scribe('python analysis/src/scatterplot_correlation_with_controls.py -i %s' % with_slash)
        # Fig 2.D: FBE vs RANK.
        os.system('echo')
        os.system('echo python analysis/src/plot_fbe_vs_rank.py -i %s' % with_slash)
        os.system('python analysis/src/plot_fbe_vs_rank.py -i %s' % with_slash)
        # Outpus to figs/method/basename_one_to_two.eps.

        # Fig 2.E Overlap with RIP-ChIP.
        os.system('echo')
        os.system('echo python analysis/src/venn_of_comparison_to_rip_chip.py %s' % with_slash)
        os.system('python analysis/src/venn_of_comparison_to_rip_chip.py %s' % with_slash)
        # Writes to figs/${method}/overlap_with_rip_*.

        # Fig 3: FBF-1 vs FBF-2.
        #     A. FBF-1/FBF-2 scatterplot correlation.
        #     B. FBF-1/FBF-2 correlation heatmap.
        #     C. patr-1 or other RNA showing differences.

        # Fig 3.A: Scatterplot of correlation between FBF-1 and FBF-2.
        # -p option here will load existing data to save time.
        scribe('python analysis/src/scatterplot_correlation_by_wig_ui.py --no_ui -d %s' % with_slash)
        scribe('python analysis/src/scatterplot_correlation_fbf1_and_fbf2.py -d %s' % with_slash)
        # Outputs to /figs/method/scatterplot_by_wig.pdf and
        # data/method/outliers_fbf1_vs_fbf2.txt.

    # Fig 3.B: Tables are output with Fig 2.2, and need GO run on them.
    os.system('echo')
    os.system('echo python analysis/src/extract_ranges_and_write_as_bed.py %s' % with_slash)
    #os.system('echo skipping...')
    os.system('python analysis/src/extract_ranges_and_write_as_bed.py %s' % with_slash)
    os.system('echo')
    #os.system('echo python analysis/src/call_bamCorrelate.py %s' % with_slash)
    #os.system('echo skipping...')
    #os.system('python analysis/src/extract_ranges_and_write_as_bed.py %s' % with_slash)
    os.system('echo')
    #os.system('python analysis/src/call_bamCorrelate.py %s' % with_slash)
    # Outputs to figs/method/spearman_correlation_between_bams_in_peak_regions.png.

    # Fig 3.C: Get from IGV.

    if is_five_reps:
        # Fig 3.4: Tables of RNA types.
        scribe('python analysis/src/table_of_rna_types.py %s' % with_slash)
        scribe('python analysis/src/boxplot_of_peaks_in_dif_locations.py -s 500 -i %s/combined_fbf.txt' % without_slash)

        scribe('python analysis/src/compare_with_ripchip.py %s' % with_slash)
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
    #os.system('echo Running CIMS and CITS analysis...')
    #os.system('echo skpping')
    #export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/
    #(cd ../raw\ fastq\ for\ Doug\ CIMS\ CITS/ && exec python ../src/cims_master.py)
    #(cd ../raw\ fastq\ for\ Doug\ CIMS\ CITS/ && exec python ../src/cims_pos_vs_fbe.py)

    # Fig 5.C: MEME-ChIP results.
    # Run on the reproducible CIMS/CITS .fasta files put in
    # raw\ fastq\ for\ Doug\ CIMS\ CITS/fasta/.

    # Fig 5.D: MEME on the peak caller output.
    # Fasta files are in data/fasta/.
    
    # Fig 5.E: Alternative FBEs.    
    #scribe('python analysis/src/alternative_fbes.py %s' % with_slash)

    # Outputs to data/method/cryptic_fbes_numofshuffles.txt.

    # Fig 6: Binding site location.
    #     A. Pie chart of binding location.
    #     B. Position in a normalized transcript.
    #     C. Peak traces for 5 candidate 5 prime UTR targets.
    #     D. Peak location in 3 prime UTRs.

    if is_five_reps:
        # Fig 6.A: Pie chart of binding site location.
        scribe('python analysis/src/peak_locations.py %s' % with_slash)
        # Outputs to figs/method/peak_locations_pie_chart_label.pdf.

        # Fig 6.B, 6.D: Position in a normalized transcript.
        scribe('python analysis/src/feature_locations/determine_feature_locations.py -i %s/combined_fbf.txt' % without_slash)
        # Outputs to figs/features_in_normalized_gene.pdf,
        # figs/features_in_utr.pdf

        # Fig 6 B, D: By subpeak.
        scribe('python analysis/src/feature_locations/subpeaks_ui.py -i {s}/combined_fbf.txt --no_ui'.format(s=without_slash))

    # Fig 6.C: Get from IGV.

    # Fig 7: Different methods of controlling iCLIP data.
    #     A. Three barplots for the methods of initial control.
    #     B. Barplot of one method vs a gradient of enrichment ratio.

    # Fig 7.A: Different methods of initial control.
    #scribe('python analysis/src/test_various_controls.py -b -g config.ini')
    #scribe('python analysis/src/test_various_controls.py --info')
    #scribe('python analysis/src/test_various_controls.py -c -g config.ini')
    # Outputs to figs/various_controls/effect_of_controls_barplot_.pdf

    #if is_five_reps:
        # Fig 7.B: Gradient of entichment vs cutoff.
        #scribe('python analysis/src/effect_of_different_controls_figure.py --input %s/combined_fbf.txt' % without_slash)
        # Outputs to figs/method/effect_of_controls_barplot_*.pdf

    # Fig SX: Bernstein figure.
    #scribe('python analysis/src/feature_locations/mutations_in_fbe.py data/fasta/subpeaks/highest_both.fa')

    # Currently unused: dCLIP.
    # dCLIP needs SAM files.
    #os.system('echo')
    #os.system('echo python analysis/src/call_dCLIP.py ${method}')
    #os.system('python analysis/src/call_dCLIP.py ${method}')
    #os.system('echo python analysis/src/get_fasta_of_dCLIP-run_DREME.py ${method}')
    #os.system('python analysis/src/get_fasta_of_dCLIP-run_DREME.py ${method}')


def scribe(cmd):
    os.system('echo')
    os.system('echo ' + cmd)
    os.system(cmd)


if __name__ == '__main__':
    for n in range(0, 7): print '**' * (3 - abs(n - 3))
    method = sys.argv[1]
    if os.path.isfile(method):
        method = os.path.dirname(method)
        print "Treating the input file as an input directory: {d}".format(
            d=method
        )
    make_filtered_sep_reps = False
    if make_filtered_sep_reps:
        os.system('mkdir for_methods_comparison/separate_reps_filtered/')
        for dirname in glob.glob('for_methods_comparison/separate_reps/*'):
            msg = ' python src/filter_peaks_by_ratio_in_peak_region.py' + \
        ' -f --analysis -i {i} --output_dir {o} --use respective '.format(
            i=dirname,
            o='for_methods_comparison/separate_reps_filtered/' + os.path.basename(dirname))
            print msg
            os.system(msg)
    make_figs(method)
    do_all = False
    if do_all:
        if not is_filtered:
            make_figs('{d}/filtered_{i}/'.format(
                d=os.path.dirname(method), i=os.path.basename(method)
            ))
        if not is_five_reps:
            make_figs('{d}/{i}_five_reps/'.format(
                d=os.path.dirname(method), i=os.path.basename(method)
            ))
        if not (is_five_reps or is_filtered):
            make_figs('{d}/filtered_{i}_five_reps/'.format(
                d=os.path.dirname(method), i=os.path.basename(method)
            ))
