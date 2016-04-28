from add_reads_columns import load_bedgraph, add_heights_to_peak_file
import pandas
import glob


def load_ga():
    ga = {}
    ga['fbf1_reads'] = load_bedgraph('bedgraph_norm/combined_fbf1.wig')
    ga['fbf2_reads'] = load_bedgraph('bedgraph_norm/combined_fbf2.wig')
    ga['n2_from_fbf1_reads'] = load_bedgraph('bedgraph_norm/combined_fbf1_n2.wig')
    ga['n2_from_fbf2_reads'] = load_bedgraph('bedgraph_norm/combined_fbf2_n2.wig')
    ga['rna_seq_reads'] = load_bedgraph(
        'bedgraph_norm/rna_seq_mod4594.wig')
    return ga


def run(indir):
    ga = load_ga()
    for fname in glob.glob(indir + '/*.txt'):
        #peaks = pandas.read_csv(fname, sep='\t')
    #    for exp in ['fbf1_reads', 'fbf2_reads', 'n2_from_fbf1_reads',
    #                'n2_from_fbf2_reads', 'rna_seq_reads']:
        add_heights_to_peak_file(fname, ga)
        #peaks.to_csv(fname, sep='\t', index=False)

if __name__== '__main__':
    import sys
    run(sys.argv[1])
