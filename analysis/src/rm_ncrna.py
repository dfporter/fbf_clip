import pandas
import os
import glob

new_dir = 'only_mrna_methods/'
os.system('mkdir ' + new_dir)
for fname in glob.glob('methods/*/*.txt'):
    print fname
    peaks = pandas.read_csv(fname, sep='\t')
    peaks = peaks[peaks['biotype']=='protein_coding']
    new_sub_dir = new_dir +'/%s' % os.path.basename(os.path.dirname(fname))
    os.system('mkdir ' + new_sub_dir)
    peaks.to_csv(new_sub_dir + '/%s' % os.path.basename(fname),
                 sep='\t', index=False)