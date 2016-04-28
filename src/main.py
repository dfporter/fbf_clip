def find_peaks(coverage, chrm='I', strand='+'):
    for iv in coverage[chrm][strand].steps():
        pass
    last_value = iv[0].start
    min_read_cutoff = 10
    peak_pos = []
    width = 1e3
    for start in range(0, last_value, int(width)):
        end = start + width
        window = HTSeq.GenomicInterval('I', start, end, strand)
        a = np.fromiter(coverage[window], dtype='i')
        if max(a) < min_read_cutoff:
            continue
        xs = np.arange(100,200,100)
        peakind = signal.find_peaks_cwt(a, xs, noise_perc=20)
        for a_pos in peakind:
            peak_pos.append(start + a_pos)
        if start > 1e6:
            break
    return peak_pos

if __name__ == '__main__':
    gtffile = HTSeq.GFF_Reader("celegans_genome/wormbase_ws235/Caenorhabditis_elegans.WBcel235.78.gtf")
    #bamfile = HTSeq.BAM_Reader( "bams/combined_fbf.bam")
    bamfile = HTSeq.BAM_Reader( "/home/dp/Desktop/bams/celegans/run813_fbf_aacc_20mapq.bam")
    coverage = HTSeq.GenomicArray("auto", stranded=True, typecode='i')
    for aln in bamfile:
        if aln.aligned:
            coverage[aln.iv] += 1
    peak_pos = find_peaks(coverage, chrm=chrm, strand=strand)
    peak_objs = []
    determine_borders(peak_pos, coverage)
    peaks_by_chrm = merge_overlapping(peak_objs)
    genes = asign_to_gene(gtffile, peaks_by_chrm)
    
