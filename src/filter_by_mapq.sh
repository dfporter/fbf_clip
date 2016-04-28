for f in ./*.bam; do
	outname=$(basename $f .bam)_20mapq.bam
	echo "samtools view -b -q 20 $f > $outname"
	#samtools view -b -q 20 $f > $outname
	sorted_out=$(basename $outname .bam)
	samtools sort $outname $sorted_out
	samtools index $outname
done
