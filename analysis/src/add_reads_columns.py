import HTSeq
import numpy as np
import pandas
import glob

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
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), strand)] = float(s[3])
    return ga


def get_val(ga, iv):
    return np.max(np.fromiter(ga[iv], dtype=np.float))


def get_pos_of_max_coverage(ga, iv):
    index = np.argmax(np.fromiter(ga[iv], dtype=np.float))
    return iv.start + index


def add_heights_to_peak_file(peak_fname, bed_ga_dict):
    try:
        peaks = pandas.read_csv(peak_fname, sep='\t')
    except:
        print "Could not open %s" % peak_fname
        return
    try:
        ivs = zip(peaks['chrm'].tolist(), peaks['left'].tolist(),
            peaks['right'].tolist(), peaks['strand'].tolist())
    except:
        print "Could not create ivs from file %s" % peak_fname
        return
    for bedgraph_name in bed_ga_dict:
        peaks[bedgraph_name] =  [
            get_val(bed_ga_dict[bedgraph_name], HTSeq.GenomicInterval(*iv)
                    ) for iv in ivs]
        if bedgraph_name in ['fbf1_reads', 'fbf2_reads']:
            peaks[bedgraph_name + '_pos_of_max_coverage'] =  [
                get_pos_of_max_coverage(bed_ga_dict[bedgraph_name], HTSeq.GenomicInterval(*iv)
                        ) for iv in ivs]
    peaks.to_csv(peak_fname, sep='\t', index=False)
    return peaks


if __name__ == '__main__':
    for fname in glob.glob('fbf1/peaks/*/*.txt'):
        print fname
    ga = {}
    ga['fbf1_reads'] = load_bedgraph('bedgraph_norm/combined_fbf1.wig')
    ga['fbf2_reads'] = load_bedgraph('bedgraph_norm/combined_fbf2.wig')
    ga['n2_from_fbf1_reads'] = load_bedgraph('bedgraph_norm/combined_fbf1_n2.wig')
    ga['n2_from_fbf2_reads'] = load_bedgraph('bedgraph_norm/combined_fbf2_n2.wig')
    ga['rna_seq_reads'] = load_bedgraph('bedgraph_norm/rna_seq_mod4594.wig')
    for fname in glob.glob('fbf1/peaks/*/*.txt'):
        print fname
        add_heights_to_peak_file(fname, ga)
    for fname in glob.glob('fbf2/peaks/*/*.txt'):
        print fname
        add_heights_to_peak_file(fname, ga)
