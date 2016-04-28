
#bowtie2 -x ../Celegans_genome/Sequence/Bowtie2Index/genome -U runx_fbf_tccg_unique_trimmed.fastq --local --phred33 -S runx_fbf_tccg.sam
#bowtie2 -x ../Celegans_genome/Sequence/Bowtie2Index/genome -U runx_fbf_aacc_unique_trimmed.fastq --local --phred33 -S runx_fbf_aacc.sam
bowtie2 -x ../Celegans_genome/Sequence/Bowtie2Index/genome -U runx_fbf_gcca_unique_trimmed.fastq --local --phred33 -S runx_fbf_gcca.sam

#bowtie2 -x ../Celegans_genome/Sequence/Bowtie2Index/genome -U runx_n2_acaa_unique_trimmed.fastq --local --phred33 -S runx_n2_acaa.sam
#bowtie2 -x ../Celegans_genome/Sequence/Bowtie2Index/genome -U runx_n2_ccgg_unique_trimmed.fastq --local --phred33 -S runx_n2_ccgg.sam
#bowtie2 -x ../Celegans_genome/Sequence/Bowtie2Index/genome -U runx_n2_tgcc_unique_trimmed.fastq --local --phred33 -S runx_n2_tgcc.sam

