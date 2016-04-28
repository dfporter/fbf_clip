import HTSeq
'''
python bed_to_wig.py input_dir output_dir

Creates bedgraph files of the read five prime position.
'''

__author__ = 'dp'
import glob
import re
import sys
import os

indir = sys.argv[1]
outdir = sys.argv[2]

ga_all_fog = HTSeq.GenomicArray('auto', stranded=True)
ga_all_control = HTSeq.GenomicArray('auto', stranded=True)

def add_to_ga(infile, global_ga):
    ga = HTSeq.GenomicArray('auto', stranded=True)
    with open(infile, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            if s[5] == '+':
                iv = HTSeq.GenomicPosition(
                    s[0], int(s[1]), s[5])
            if s[5] == '-':
                iv = HTSeq.GenomicPosition(
                    s[0], int(s[1]) - 1, s[5])
            ga[iv] += 1
            global_ga[iv] += 1
    return ga

def add_from_wig(infile, global_ga):
    if re.match('.*_-_.*', os.path.basename(infile)) is not None:
        strand = '-'
    if re.match('.*_\+_.*', os.path.basename(infile)) is not None:
        strand = '+'
    with open(infile, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            iv = HTSeq.GenomicInterval(
                    s[0], int(s[1]), int(s[2]), strand)
            global_ga[iv] += float(s[3])
    return global_ga

ga_all_fog_coverage = HTSeq.GenomicArray('auto', stranded=True)
ga_all_control_coverage = HTSeq.GenomicArray('auto', stranded=True)

for infile in glob.glob(indir +'*.wig'):
    if re.match('.*fog.*', os.path.basename(infile)) is not None:
        add_from_wig(infile, ga_all_fog_coverage)
    if re.match('.*control.*', os.path.basename(infile)) is not None:
        add_from_wig(infile, ga_all_control_coverage)

ga_all_fog_coverage.write_bedgraph_file('wig/all_fog_+.wig', strand='+')
ga_all_fog_coverage.write_bedgraph_file('wig/all_fog_-.wig', strand='-')

sys.exit()
for infile in glob.glob(indir +'*.bed'):
    if re.match('.*fog.*', os.path.basename(infile)) is not None:
        ga = add_to_ga(infile, ga_all_fog)
    if re.match('.*control.*', os.path.basename(infile)) is not None:
        ga = add_to_ga(infile, ga_all_control)
    outname = "{d}/{b}".format(
        d=outdir,
        b=os.path.basename(infile).partition('.bed')[0])
    outname_plus = outname + '_+.wig'
    ga.write_bedgraph_file(outname_plus, strand='+')
    outname_minus = outname + '_-.wig'
    ga.write_bedgraph_file(outname_minus, strand='-')

ga_all_fog.write_bedgraph_file(outdir + '/all_fog_+.wig', strand='+')
ga_all_fog.write_bedgraph_file(outdir + '/all_fog_-.wig', strand='-')
ga_all_control.write_bedgraph_file(outdir + '/all_control_+.wig', strand='+')
ga_all_control.write_bedgraph_file(outdir + '/all_control_-.wig', strand='-')
