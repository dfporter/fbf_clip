"""

For a directory of peaks files, e.g. pypeaks_fdr1_negip/,
given with -i, filter every peaks file by the ratio of reads in the
peak region over a control, and output to a directory given with -o.

"""

import sys
sys.path.insert(0, '/Network/Servers/file.biochem.wisc.edu/Volumes/BioSAN/Users/dfporter/.local/lib/python2.6/site-packages')
sys.path.insert(0, './src/')
sys.path.insert(0, './analysis/src/')
import HTSeq
import pandas
import os
import argparse
import subprocess
import callpeaks
import glob
from effect_of_different_controls_figure import *
from filter_peaks_by_ratio_in_peak_region import *

def load_bed_info_for_a_directory(dir_name):
    """Overwrites input files with versions with peak num columns.
    """
    print "load_bed_info called with dir_name %s" % dir_name
    ga = {}
    ga_five = {}
    (ga['fbf1_reads'], ga_five['fbf1_reads']) = get_bed(
        'fbf1/beds/combined_fbf1.bed')
    (ga['fbf2_reads'], ga_five['fbf2_reads']) = get_bed(
        'fbf2/beds/combined_fbf2.bed')
    (ga['n2_from_fbf1_reads'], ga_five['n2_from_fbf1_reads']) = get_bed(
        'fbf1/beds/runx_n2_combined_20mapq.bed')
    (ga['n2_from_fbf2_reads'], ga_five['n2_from_fbf2_reads']) = get_bed(
        'fbf2/beds/run813_n2_combined_20mapq.bed')
    (ga['rna_seq_reads'], ga_five['rna_seq_reads']) = get_bed(
        'lib/rna_seq/modencode_4594.bed')
    for filename in glob.glob(dir_name + '/*.txt'):
        peaks = pandas.DataFrame.from_csv(
            filename, sep='\t', index_col=False)
        print "Loaded peaks from filename %s" % filename
        print "Would output to rewrite the same file."
        for label in ga_five.keys():
            add_reads_in_peak(peaks, ga_five[label], label=label)
        peaks.to_csv(filename, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__)
    parser.add_argument(
        '-d', '--dont_overwrite',
        dest='dont_overwrite',
        default=True,
        help="""If a given column label already exists, don't overwrite it.""")
    parser.add_argument('-i', '--input', dest='input',
                        default='with_peak_nums/combined_fbf_4_rna_seq.txt',
                        help='Input directory (Required).')
    parser.add_argument('-o', '--output',
                        default='filtered_peaks/',
                        help="Directory to output filtered peaks to.")
    args = parser.parse_args()
    if not os.path.exists(args.output + '/'):
        os.system('mkdir ' + args.output)
    for peaks_filename in glob.glob(args.input + '/*.txt'):
        print 'Filtering %s...' % peaks_filename
        out_filename = args.output
        out_filename += '/' + os.path.basename(peaks_filename)
        print 'Outputting to %s...' % out_filename
        peaks = pandas.DataFrame.from_csv(
            peaks_filename, sep='\t', index_col=False)
        need_info = False
        for colname in ['fbf1_reads', 'n2_from_fbf1_reads',
                        'fbf2_reads', 'n2_from_fbf2_reads']:
            if colname not in peaks.columns:
                need_info = True
                print "need info col %s for %s" % (colname, peaks_filename)
        #if need_info:
        #    load_bed_info_for_a_directory(args.input)
        #sys.exit()
        peaks = convert_df_to_list_of_dicts(peaks)
        filtered_peaks = filter_list_by_ratio(
            peaks, 'fbf1_reads', 'n2_from_fbf1_reads', 10)
        filtered_peaks = filter_list_by_ratio(
            filtered_peaks, 'fbf2_reads', 'n2_from_fbf2_reads', 10)
        filtered_peaks_df = pandas.DataFrame(filtered_peaks)
        filtered_peaks_df.to_csv(
            out_filename, index=False, sep='\t')

