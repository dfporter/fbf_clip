import pandas
import re


def make_gtf_with_biotype_column(
    gtf_filename='../clip/lib/gtf_with_names_column.txt',
    out_filename='../clip/lib/gtf_with_biotype_column2.txt',
    only_cDNA=False):
    gtf_df = pandas.read_csv(gtf_filename, sep='\t')
    gtf_df['biotype'] = ''
    for index, row in gtf_df.iterrows():
        if index % 1000 == 0: print index
        m = re.search('gene_biotype "([^"]+)"', gtf_df.loc[index, '8'])
        if m is not None:
            gtf_df.loc[index, 'biotype'] = m.group(1)
    if only_cDNA:
        gtf_df = gtf_df[gtf_df['biotype']=='protein_coding']
    gtf_df.to_csv(out_filename, sep='\t', index=False)


if __name__ == '__main__':
    make_gtf_with_biotype_column()

