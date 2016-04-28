import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
import matplotlib
import os
import xlwt


def categorize_rna_types(peaks, label='', subdir=''):
    biotypes = {}
    if 'biotype' not in peaks.columns:
        for index, row in peaks.iterrows():
            m = re.search('gene_biotype "([^"]+)"', row['8'])
            if m is not None:
                biotypes.setdefault(m.group(1), set())
                biotypes[m.group(1)].add(row['gene_name'])
    else:
        for index, row in peaks.iterrows():
            biotypes.setdefault(row['biotype'], set())
            biotypes[row['biotype']].add(row['gene_name'])
    li = "Experiment: %s\n" % label
    for _type in biotypes:
        li += "\t%s: %i" % (_type, len(list(biotypes[_type])))
        if len(list(biotypes[_type])) <= 100:
            li += "\t" + ", ".join(biotypes[_type])
        li += "\n"
    print li
    print "***"
    book = xlwt.Workbook()
    sh = book.add_sheet('stats')
    sh.write(0, 0, 'Gene type')
    sh.write(0, 1, 'No. targets')
    sh.write(0, 2, 'Examples')
    for n, _type in enumerate(biotypes, start=1):
        sh.write(n, 0, _type)
        sh.write(
            n, 1, len(list(biotypes[_type]))
            )
        if len(list(biotypes[_type])) <= 100:
            examples = ", ".join(biotypes[_type])
        else: examples = ""
        sh.write(n, 2, examples)
    if not os.path.exists('tables/'):
        os.system('mkdir tables/')
    if not os.path.exists('tables/%s/' % subdir):
        os.system('mkdir ' + 'tables/%s/' % subdir)
    book.save('./tables/%s/rna_types_table_%s.xls' % (subdir, label))
    
if __name__ == '__main__':
    gtf_sep_cols = pandas.read_csv('lib/gtf_with_names_column.txt', sep='\t')
    #top_level_dir = './pypeaks_fdr5_negip_local/'
    top_level_dir = sys.argv[1]
    peaks = {}
    for filename in glob.glob(top_level_dir + '/combined_*.txt'):
        peaks[filename] = pandas.read_csv(filename, sep='\t')
    for exp in peaks:
        #pie_chart_of_peak_locations(peaks[exp], label=os.path.basename(exp),
        #                            subdir=os.path.dirname(top_level_dir))
        peaks[exp].sort('height', inplace=True, ascending=False)
        peaks[exp].drop_duplicates('gene_name', inplace=True)
    for exp in peaks:
        categorize_rna_types(peaks[exp], label='highest_peak_per_gene_%s' % os.path.basename(exp),
                                    subdir=os.path.dirname(top_level_dir))        
