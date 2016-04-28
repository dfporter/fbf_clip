import pandas
import re

def map_type(name, gene_types):
    print "Looking for -%s-" % name 
    if name in gene_types:
        print "Found: %s" % gene_types[name]
        return gene_types[name]
    else: return "Unknown"

def annotate_peaks_with_gene_type(
        peaks,
        gtf_filename='../clip/lib/gtf_with_biotype_column.txt',
        output_filename=False,
        only_cDNA=False, gtf_df=False):
    print "Called annotate_peaks_with_gene_type"
    print peaks.head()
    if not gtf_df:
        gtf_df = pandas.read_csv(gtf_filename, sep='\t')
        no_dup = gtf_df.drop_duplicates('gene_name', inplace=False)
    gene_types = dict(
        zip(no_dup['gene_name'].tolist(), no_dup['biotype'].tolist())
    )
    peaks['biotype'] = ''
    if len(peaks.index) == 0: return
    peaks['biotype'] = peaks.apply(lambda row: map_type(row['gene_name'], gene_types),
                axis=1)


def remove_ncRNA(peaks):
    peaks = peaks[peaks['_biotype']=='protein_coding']

if __name__ == "__main__":
    peaks = pandas.read_csv('cims_tables/combined_fbf1.txt', sep='\t')
    annotate_peaks_with_gene_type(peaks)
