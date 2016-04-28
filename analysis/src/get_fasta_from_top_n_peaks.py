import glob
import sys
import re
import pandas
import os
fname = sys.argv[1]

peaks = pandas.read_csv(fname, sep='\t')
peaks.sort('height', inplace=True, ascending=False)
outli = ''
seq_dict = {}
w_fbe = ''
wo_fbe = ''
rank_cutoff = 2000
rank = 0
peaks_l = peaks.to_dict('records')
for index, row in enumerate(peaks_l):
    if re.match('.*' + str('[aA]') * 10 + '.*', row['seq']) is not None:
        if len(row['seq']) < 50:
            print 'skipping ' + row['seq']
            continue
    rank += 1
    if rank > rank_cutoff:
        break
    label = '>{g}_{n}'.format(g=row['gene_name'], n=index)
    seq_dict[label] = row['seq']

for label in seq_dict:
    seq = seq_dict[label]
    if re.match('.*[tT][gG][tT]\w\w\w[Aa][Tt].*', seq) is not None:
        w_fbe += label + '\n' + seq + '\n'
    else:
        wo_fbe += label + '\n' + seq + '\n'
print "Number of peaks:"
print len(seq_dict)
print "Top peaks:"
print peaks.head()
with open('data/fasta/subpeaks/{b}_{t}_w_fbe.fa'.format(
        b=os.path.basename(fname).partition('.txt')[0],
        t=rank_cutoff), 'w') as f:
    f.write(w_fbe)
with open('data/fasta/subpeaks/{b}_{t}_wo_fbe.fa'.format(
        b=os.path.basename(fname).partition('.txt')[0],
        t=rank_cutoff), 'w') as f:
    f.write(wo_fbe)