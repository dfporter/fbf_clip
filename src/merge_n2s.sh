#samtools merge -h runx_n2_acaa_20mapq.bam runx_n2_combined_20mapq.bam runx_n2_acaa_20mapq.bam runx_n2_tgcc_20mapq.bam runx_n2_ccgg_20mapq.bam

samtools sort runx_n2_combined_20mapq.bam runx_n2_combined_20mapq
samtools index runx_n2_combined_20mapq.bam