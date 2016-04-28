import glob
import re
import sys
import os
import subprocess


def get_total_reads(bam_filename):
        cmdl = ['samtools', 'flagstat', bam_filename]
        statsS = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
        statsS.wait()
        stats_out = statsS.communicate()[0].rstrip('\n')
        for li in stats_out.split('\n'):
                m = re.match(r'(\d+) \+ \d* mapped \(([0-9\.]+).*', li)
                if(m is not None):
                        return int(m.group(1))
        return False

for bamfile in glob.glob('fbf*/bams/*.bam'):
        #if (re.search('parclip', bamfile) is None) and (
        #        re.search('combined_clip', bamfile) is None):
        #        continue
        # Get the total dataset size
        total_reads = get_total_reads(bamfile)
        if(not total_reads):
                print "Error in getting total read number for %s" % bamfile
                continue
        scaling_factor = float('1e6')/float(total_reads)
        #scaling_factor = 1.0
        print "%i total reads in %s... Scaling by %f" % (total_reads,
                                                         bamfile, scaling_factor)
        bam_basename = os.path.basename(bamfile)
	wig_dir = re.sub('bams', 'wigs', os.path.dirname(bamfile))
	print 'Wig directory (created if not existant): ' + wig_dir
	if not os.path.exists(wig_dir): os.system('mkdir ' + wig_dir)
	for strand in ['+', '-']:
		wigfile = wig_dir +'/' + bam_basename.partition('.bam')[0] + '_%s.wig' % strand
		cmd = "genomeCoverageBed -ibam %s -bg -strand %s -scale %s" % (bamfile, strand, scaling_factor)
		cmd += " -g ../ensemblEF4/chr_lengths > %s""" % (wigfile)
		print cmd
		os.system(cmd)

