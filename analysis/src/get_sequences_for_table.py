import pandas
import HTSeq

def complement(s):
    if s[0].isupper():
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    if s[0].islower():
        print "%s is lower case." % s[0]
        basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s

def get_sequences_for_table(peaks, sequences, expand=30, max_width=False,
                            asymmetric_expand=False, left_expand=False,
                            right_expand=False):
    peaks['seq'] = ''
    peaks.left.apply(int)
    peaks.right.apply(int)
    for index in peaks.index:
        if asymmetric_expand:
            seq = sequences[peaks.loc[index, 'chrm']][
                peaks.loc[index, 'left']-left_expand:peaks.loc[index, 'right']+right_expand]
        else:
            try:
                left_bound = int(peaks.loc[index, 'left']-expand)
                right_bound = int(peaks.loc[index, 'right']+expand)
            except:
                print "Could not cast boundaries to int: %s" % str(peaks.loc[index])
            seq = sequences[peaks.loc[index, 'chrm']][
                left_bound:right_bound]
        if peaks.loc[index, 'strand'] == '-':
            seq = rc(seq)
        peaks.loc[index, 'seq'] = seq
    peaks['peak_len'] = 0
    for index in peaks.index:
        peaks.loc[index, 'peak_len'] = len(peaks.loc[index, 'seq'])
    peaks = peaks[peaks['peak_len'] > 9]
