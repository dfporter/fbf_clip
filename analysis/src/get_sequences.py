import HTSeq
import re
import pandas
from rc import rc
def get_sequences(combined, fasta_filename=None):
    """combined = pandas.DataFrame."""
    if fasta_filename is None:
        fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    if 'chrm' in combined.columns: chrm_key = 'chrm'
    elif 'chrom' in combined.columns: chrm_key = 'chrom'
    elif 'Chr' in combined.columns: chrm_key = 'Chr'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for index, peak_row in combined.iterrows():
        start = combined.loc[index, 'left']
        end = combined.loc[index, 'right']
        chrm = combined.loc[index, chrm_key]
        seq = sequences[chrm][start:end]
        if combined.loc[index, 'strand'] == '-':
            seq = rc(seq)
        combined.loc[index, 'seq'] = seq
