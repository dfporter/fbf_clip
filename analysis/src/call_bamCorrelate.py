import sys
import os
import re
import glob

bed_filename = sys.argv[1]
bam_files = [
    'fbf1/bams/runx_fbf_aacc_20mapq.bam',
    'fbf1/bams/runx_fbf_gcca_20mapq.bam',
    'fbf1/bams/runx_fbf_tccg_20mapq.bam',
    'fbf1/bams/runx_n2_acaa_20mapq.bam',
    'fbf1/bams/runx_n2_ccgg_20mapq.bam',
    'fbf1/bams/runx_n2_tgcc_20mapq.bam',
    'fbf2/bams/run813_fbf_aacc_20mapq.bam',
    'fbf2/bams/run813_fbf_gcca_20mapq.bam',
    'fbf2/bams/run813_fbf_tccg_20mapq.bam',
    'fbf2/bams/run813_n2_acaa_20mapq.bam',
    'fbf2/bams/run813_n2_ccgg_20mapq.bam',
    'fbf2/bams/run813_n2_tgcc_20mapq.bam']

bam_files_only_fbf = [
    'fbf1/bams/runx_fbf_aacc_20mapq.bam',
    'fbf1/bams/runx_fbf_gcca_20mapq.bam',
    'fbf1/bams/runx_fbf_tccg_20mapq.bam',
    'fbf2/bams/run813_fbf_aacc_20mapq.bam',
    'fbf2/bams/run813_fbf_gcca_20mapq.bam',
    'fbf2/bams/run813_fbf_tccg_20mapq.bam']

cmd = 'bamCorrelate BED-file --BED {ib} -b '.format(ib=bed_filename)

for bam_file in bam_files:
    cmd += ' %s' % bam_file

figdir = 'figs/' + os.path.basename(bed_filename.partition('.bed')[0])  + '/'
if not os.path.exists(figdir):
    os.system('mkdir %s' % figdir)
cmd += ' -o {o} --corMethod spearman'.format(
    o=figdir + 'spearman_correlation_between_bams_in_peak_regions.png')

print cmd
os.system(cmd)


cmd = 'bamCorrelate BED-file --BED {ib} -b '.format(ib=bed_filename)

for bam_file in bam_files_only_fbf:
    cmd += ' %s' % bam_file

figdir = 'figs/' + os.path.basename(bed_filename.partition('.bed')[0])  + '/'
if not os.path.exists(figdir):
    os.system('mkdir %s' % figdir)
cmd += ' -o {o} --corMethod spearman'.format(
    o=figdir + 'spearman_correlation_between_bams_in_peak_regions_only_fbf.png')

print cmd
os.system(cmd)
