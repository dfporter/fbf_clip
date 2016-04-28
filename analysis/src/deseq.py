import pandas
import os
import sys
import collections

fname = 'combined_counts.txt'
l1r = ''
l2r = ''
outli = ''
with open(fname) as f:
    l1r = next(f)
    l2r = next(f)
    for li in f:
        outli += li
if len(l1r.split('\t')) == (len(l2r.split('\t')) -1):
    l1r = 'gene_name\t' + l1r
    with open('counts/fix_header_combined_counts.txt', 'w') as f:
        f.write(l1r)
        f.write(l2r)
        f.write(outli)
    print "Fixed header to %s" % l1r
    fname = 'counts/fix_header_combined_counts.txt'
    print "Now opening the fixed file %s" % fname

cf = pandas.read_csv(, sep='\t')


rscript = """
library(DESeq2)
get_res = function(directory) {
        sf = grep('txt', list.files(directory), value=TRUE)
        sampleCondition = sub("(.*)_.*", "\\1", sf)
        sampleTable = data.frame(sampleName=sf, fileName=sf, condition=sampleCondition)
        ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,  directory=directory, design= ~ condition)
        ddsHTSeq = DESeq(ddsHTSeq)
        dds = ddsHTSeq
        res = results(dds)
        resO = res[order(res$padj),]
        return(resO)
}
in_directory = "counts/"


res_wt = get_res('counts/')
write.table(res_wt, file='enriched.txt', quote=FALSE, sep='\t')
"""
