from peaks import peaks
import HTSeq
import numpy as np
import os

def load_bedgraph(fname):
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    plus_file = fname.partition('.wig')[0] + '_+.wig'
    add_strand_to_ga_from_bedgraph_file(plus_file, ga, '+')
    minus_file = fname.partition('.wig')[0] + '_-.wig'
    add_strand_to_ga_from_bedgraph_file(minus_file, ga, '-')
    return ga


def add_strand_to_ga_from_bedgraph_file(fname, ga, strand):
    with open(fname, 'r') as f:
        next(f)
        for line in f:
            #if line[0] != 'M': continue
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), strand)] = float(s[3])
    return ga


def get_val(ga, iv):
    return np.max(np.fromiter(ga[iv], dtype=np.float))


def get_pos_of_max_coverage(ga, iv):
    index = np.argmax(np.fromiter(ga[iv], dtype=np.float))
    return iv.start + index


def run():
    rep = {}
    ga = load_bedgraph('fbf1/bedgraph_norm/control_n2.wig')
    # Unfiltered.
    infile = 'unfiltered/combined_fbf1.txt'
    fbf1 = peaks(file=infile, name='fbf1')
    fbf1.add_reads(ga=ga, name='fbf1_n2')
    fbf1.data['n2_from_fbf1_reads'] = fbf1.data['fbf1_n2']
    fbf1.data['fbf1_n2_reads'] = fbf1.data['fbf1_n2']
    fbf1.write_table(infile)
    # Filtered.
    infile = 'filtered/combined_fbf1.txt'
    fbf1 = peaks(file=infile, name='fbf1')
    fbf1.add_reads(ga=ga, name='fbf1_n2')
    fbf1.data['n2_from_fbf1_reads'] = fbf1.data['fbf1_n2']
    fbf1.data['fbf1_n2_reads'] = fbf1.data['fbf1_n2']
    fbf1.write_table(infile)
    for subdir, dirs, files in os.walk('fbf1_to_fbf2_n2/peaks/'):
        for file in files:
            path = os.path.join(subdir, file)
            p = peaks(file=path)
            p.add_reads(ga=ga, name='fbf1_n2')
            p.data['n2_from_fbf1_reads'] = fbf1.data['fbf1_n2']
            p.data['fbf1_n2_reads'] = fbf1.data['fbf1_n2']
            p.write_table(path)
    
if __name__ == '__main__':
    run()
