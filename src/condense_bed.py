import os
import sys
import re
import glob

indir = sys.argv[1] + '/'
os.system('ls ' + indir)

outdir = 'compressed_bed/'
if not os.path.exists(outdir):
	os.system('mkdir '+ outdir)

def read_bed(fname):
	with open(fname, 'r') as f:
		outli = ""
		for read_n, li in enumerate(f):
			s = li.rstrip('\n').split('\t')
			outli += '\t'.join([
				s[0], s[1], s[2], str(read_n), s[4], s[5]]) + '\n'
	return outli

for fname in glob.glob(indir + '*.bed'):
	print "Reading " + fname
	outli = read_bed(fname)
	outfname = outdir + os.path.basename(fname)
	with open(outfname, 'w') as f:
		f.write(outli)


