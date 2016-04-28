import pandas
import glob

def add_biotype_to_all(
    files,
    gtf='/groups/Kimble/Aman Prasad/clip2/lib/gtf_with_names_column.txt'):
    gtf = pandas.read_csv(gtf, sep='\t')
    to_biotype = {} #collections.defaultdict(set)
    for tup in zip(gtf.gene_name, gtf.biotype):
        if (tup[0] in to_biotype) and (
            to_biotype[tup[0]] == 'protein_coding'):
            continue
        to_biotype[tup[0]] = tup[1]
    for fname in files:
        peaks = pandas.read_csv(fname, sep='\t')
        peaks['biotype'] = [to_biotype for x in peaks.gene_name]
        peaks.to_csv(fname, sep='\t', index=False)

