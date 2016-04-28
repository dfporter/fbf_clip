for f in ./*.sam; do
	echo ${f}
	bamname=$(basename "${f}" .sam).bam
	#echo "$bamname"
	echo "samtools view -Sb ${f} > ${bamname}"
	samtools view -Sb ${f} > ${bamname}
	bambasename=$(basename "${bamname}" .bam)
	samtools sort ${bamname} $bambasename
	samtools index ${bamname}
done
