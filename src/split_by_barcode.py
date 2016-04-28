import re
import sys
import HTSeq


def split_by_barcode(sample1_filename,
                     sample2_filename,
                     sample3_filename,
                     missing_barcode_filename,
                     initial_filename):
    sample1f = open(sample1_filename, 'w')
    sample2f = open(sample2_filename, 'w')
    sample3f = open(sample3_filename, 'w')
    missingf = open(missing_barcode_filename, 'w')
    fastq_file = HTSeq.FastqReader(initial_filename)
    total_reads = 0
    sample1_num = 0
    sample2_num = 0
    sample3_num = 0
    found = 0
    for read in fastq_file:
        total_reads += 1
	if(not (total_reads % 100000)):
		sum_with_barcode = sample1_num + sample2_num +sample3_num
		print "Read: %i (%i, %i, %i). Barcode in %f reads" % (
			total_reads,
			sample1_num,
			sample2_num,
			sample3_num,
			float(sum_with_barcode)/float(total_reads))
	if strain == 'n2':
            if(re.match('\w{3}TTGT.*', read.seq)):
                sample1_num += 1
                read.write_to_fastq_file(sample1f)
                continue
            if(re.match('\w{3}CCGG.*', read.seq)):
                sample2_num += 1
                read.write_to_fastq_file(sample2f)
                continue
            if(re.match('\w{3}GGCA.*', read.seq)):
                sample3_num += 1
                read.write_to_fastq_file(sample3f)
                continue
            read.write_to_fastq_file(missingf)
        if strain == 'fbf':
            if(re.match('\w{3}TGGC.*', read.seq)):
                sample1_num += 1
                read.write_to_fastq_file(sample1f)
                continue
            if(re.match('\w{3}GGTT.*', read.seq)):
                sample2_num += 1
                read.write_to_fastq_file(sample2f)
                continue
            if(re.match('\w{3}CGGA.*', read.seq)):
                sample3_num += 1
                read.write_to_fastq_file(sample3f)
                continue
            read.write_to_fastq_file(missingf)
    print """Results:
Reads: %i
Sample 1: %i (%f)
Sample 2: %i (%f)
(sum sample 1 + 2): %i
remaining: %i""" % (total_reads,
                    sample1_num, float(sample1_num)/float(total_reads),
                    sample2_num, float(sample2_num)/float(total_reads),
                    sample1_num + sample2_num,
                    total_reads - sample1_num + sample2_num)

if __name__ == '__main__':
    initial_filename = sys.argv[1]
    strain = sys.argv[2] # 'fbf'
    run_id = sys.argv[3] #'runx'
    split_by_barcode(
            '%s_%s_barcode1.fastq' % (run_id, strain),
            '%s_%s_barcode2.fastq' % (run_id, strain),
            '%s_%s_barcode3.fastq' % (run_id, strain),
            '%s_%s_unknown.fastq' % (run_id, strain),
            initial_filename)
    """
cgga > tccg
ggtt > aacc
tggc > accg
ggca > tgcc -> fbf2 has this barcode, fbf1 does not.

n2
ccgg > y
ttgt > fbf1 huge
ttaa > nothing for fbf1 or fbf2 n2s, really
aata > nothing for either


"""
