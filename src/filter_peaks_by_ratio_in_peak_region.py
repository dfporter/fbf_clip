import sys
sys.path.insert(
    0, '/Network/Servers/file.biochem.wisc.edu/Volumes/BioSAN/Users/dfporter/.local/lib/python2.6/site-packages')
import HTSeq
import pandas
import os
import argparse
import subprocess
import callpeaks
import glob
import re
import time
import _filter

def get_bed(bed_filename, debug_mode=False):
    if debug_mode:
        os.system('head -n 10000 %s > %s' % (bed_filename, 'tmp.bed'))
        f = open('tmp.bed', 'r')
    else:
        f = open(bed_filename, 'r')
    #ga = HTSeq.GenomicArray(chroms='auto')
    ga_five = HTSeq.GenomicArray(chroms='auto')
    for n, read in enumerate(f):
        if not n % 1e6: print "get_bed({fn}) read {n}".format(fn=bed_filename, n=n)
        s = read.rstrip('\n').split('\t')
        #iv = HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), s[5])
        #ga[iv] += 1
        if s[5] == '-':
            # Check if there should be a 1 offset here.
            ga_five[HTSeq.GenomicPosition(s[0], int(s[2]), s[5])] += 1
        if s[5] == '+':
            ga_five[HTSeq.GenomicPosition(s[0], int(s[1]), s[5])] += 1
    f.close()
    return ga_five


def add_reads_in_peak(
        peaks, ga_five, dont_overwrite=False, label='fbf_reads'):
    if dont_overwrite and (label in peaks.columns):
        return
    peaks = _filter.add_heights_to_peak_file(peaks, ga_five)
    return peaks
    peaks[label] = 0
    peaks_l = peaks.to_dict('records')
    for row in peaks_l:
        peak_iv = HTSeq.GenomicInterval(
            row['chrm'],
            row['left'],
            row['right'],
            row['strand'])
        for iv, score in ga_five[peak_iv].steps():
            row[label] += score
    peaks = pandas.DataFrame(peaks_l)
    return peaks

def add_reads_peak_using_samtools(peaks, bam_filename, label):
    for index, row in peaks.iterrows():
        cmdl = ['samtools', 'view', str(bam_filename),
                row['chrm'], ':', str(row['left']),
                '-', str(row['right'])]
        samtools_out = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
        print samtools_out.communicate()


def filter_by_ratio(peaks, col1, col2, min_1_to_2_ratio):
    assert type(peaks) is type([])
    if len(peaks) > 0: assert type(peaks[0]) is type({})
    print "\t\tFiltering %s vs %s by %f" % (col1, col2, min_1_to_2_ratio)
    for row in peaks:
        ratio = float(row[col1])/max([float(row[col2]), 1.])
        if ratio >= min_1_to_2_ratio:
            row['keep'] = 1
            if row[col1] < 10.:
                print "Set %s wrongly" % str(row)
        else:
            row['keep'] = 0
    peaks = filter(lambda row: row['keep'], peaks)
#    peaks = peaks[peaks['keep'] == 1]
    return peaks


def load_bed_and_add_to_peaks(peaks, args, _bed, label, out_filename):
    if (label not in peaks.columns) or (not args.dont_overwrite):
        (ga, ga_five) = get_bed(_bed)
        peaks = add_reads_in_peak(peaks, ga_five, label=label)
        peaks.to_csv(out_filename, sep='\t', index=False)


def add_bed_to_peaks(peaks, args, ga_five, label, out_filename):
    peaks = add_reads_in_peak(peaks, ga_five, label=label)


def do_various_filters(
        peaks, peaks_filename, out_filename, args, ratio=10., only_filter_by='both'):
    # Remove any duplicate indexes.
    #peaks['index_n'] = range(0, len(peaks.index))
    #peaks.set_index('index_n')
    print only_filter_by
    if type(peaks) is pandas.DataFrame:
        peaks = peaks.to_dict('records')
    print "\t%s: initial len %i" % (out_filename, len(peaks))
    if only_filter_by == 'both' or only_filter_by == 'fbf1':
        peaks = filter_by_ratio(peaks, 'fbf1_reads', 'n2_from_fbf1_reads', ratio)
        print "\tAfter FBF-1 filtering %i" % len(peaks)
    if only_filter_by == 'both' or only_filter_by == 'fbf2':
        peaks = filter_by_ratio(peaks, 'fbf2_reads', 'n2_from_fbf2_reads', ratio)
        print "\tAfter FBF-2 filtering %i" % len(peaks)
    if only_filter_by == 'fbf2_n2':
        if args.filter_by_both_fbfs:
            peaks = filter_by_ratio(peaks, 'fbf1_reads', 'n2_from_fbf2_reads', ratio)
        peaks = filter_by_ratio(peaks, 'fbf2_reads', 'n2_from_fbf2_reads', ratio)
    if only_filter_by == 'fbf1_n2':
        peaks = filter_by_ratio(peaks, 'fbf1_reads', 'n2_from_fbf1_reads', ratio)
        if args.filter_by_both_fbfs:
            peaks = filter_by_ratio(peaks, 'fbf2_reads', 'n2_from_fbf1_reads', ratio)
    peaks = pandas.DataFrame(peaks, index=None)
    peaks.to_csv(out_filename, sep='\t', index=False)


def load_all_bedgraph_info(args, config):
    import _filter
    ga = {}
    ga['fbf1_reads'] = _filter.load_bedgraph_file(
        'bedgraph_unnorm/combined_fbf1')
    ga['fbf2_reads'] = _filter.load_bedgraph_file(
        'bedgraph_unnorm/combined_fbf2')
    ga['n2_from_fbf1_reads'] = _filter.load_bedgraph_file('bedgraph_unnorm/fbf1_n2')
    ga['n2_from_fbf2_reads'] = _filter.load_bedgraph_file('bedgraph_unnorm/fbf2_n2')
    return ga


def load_all_bed_info(only_use='fbf1', configfile=None):
    if not os.path.exists('tmp_peaks_with_read_nums/'):
        os.system('mkdir tmp_peaks_with_read_nums/')
    ga_five = {}
    if configfile is not None:
        pass  # Not implemented.
    else:
        fbf1_fname = 'fbf1/beds/combined_fbf1.bed'
        fbf2_fname = 'fbf2/beds/combined_fbf2.bed'
        fbf1_five = 'fbf1/beds/runx_n2_combined_20mapq.bed'
        fbf2_five = 'fbf2/beds/run813_n2_combined_20mapq.bed'
        rnaseq = 'lib/rna_seq/modencode_4594.bed'
    if only_use == 'both' or only_use == 'fbf1':
        ga_five['fbf1_reads'] = get_bed(
                'fbf1/beds/combined_fbf1.bed')
    if only_use == 'both' or only_use == 'fbf2':
        ga_five['fbf2_reads'] = get_bed(
            'fbf2/beds/combined_fbf2.bed')
    if only_use == 'both' or only_use == 'fbf1':
        ga_five['n2_from_fbf1_reads'] = get_bed(
            'fbf1/beds/runx_n2_combined_20mapq.bed')
    if only_use == 'both' or only_use == 'fbf2':
        ga_five['n2_from_fbf2_reads'] = get_bed(
            'fbf2/beds/run813_n2_combined_20mapq.bed')
    ga_five['rna_seq_reads'] = get_bed(
        'lib/rna_seq/modencode_4594.bed')
    return ga_five
    # While we're here, save the read locations for later.
    #out_dir = 'data/wigs_five_prime/'
    #if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
    #for key in ga_five:
    #    ga_five[key].write_bedgraph_file(out_dir + key + '_plus.bed', '+')
    #    ga_five[key].write_bedgraph_file(out_dir + key + '_minus.bed', '-')
    #out_dir = 'data/wigs_coverage/'
    #if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
    #for key in ga:
    #    ga[key].write_bedgraph_file(out_dir + key + '_plus.bed', '+')
    #    ga[key].write_bedgraph_file(out_dir + key + '_minus.bed', '-')


def filter_using_argument(peaks, peaks_filename, out_filename, args):
    peaks = peaks.to_dict('records')
    peaks_basename = os.path.basename(peaks_filename)
    do_various_filters(
            peaks, peaks_basename, out_filename, args,
            ratio=args.ratio, only_filter_by=args.use)
    peaks = pandas.DataFrame(peaks, index=None)
    return peaks


def filter_using_filename(peaks, peaks_filename, out_filename, args):
    peaks_basename = os.path.basename(peaks_filename)
    if args.use == 'respective':
        if type(peaks) is pandas.DataFrame:
            peaks = peaks.to_dict('records')
        if re.search('fbf1', peaks_basename) is not None:
            peaks = filter_by_ratio(peaks, 'fbf1_reads', 'n2_from_fbf1_reads', args.ratio)
            print "\t\tLength %i" % len(peaks)
        if re.search('fbf2', peaks_basename) is not None:
            peaks = filter_by_ratio(peaks, 'fbf2_reads', 'n2_from_fbf2_reads', args.ratio)
            print "\t\tLength %i" % len(peaks)
        if peaks_basename == 'combined_fbf.txt':
            peaks = filter_by_ratio(peaks, 'fbf1_reads', 'n2_from_fbf1_reads', args.ratio)
            peaks = filter_by_ratio(peaks, 'fbf2_reads', 'n2_from_fbf2_reads', args.ratio)
        peaks = pandas.DataFrame(peaks, index=None)
        peaks.to_csv(out_filename, sep='\t', index=False)
        return
    if re.search('fbf1', peaks_basename) is not None:
        do_various_filters(
            peaks, peaks_basename, out_filename, args,
            ratio=args.ratio, only_filter_by='fbf1_n2')
    elif re.search('fbf2', peaks_basename) is not None:
        do_various_filters(
            peaks, peaks_basename, out_filename, args,
            ratio=args.ratio, only_filter_by='fbf2_n2')
    else:
        do_various_filters(
            peaks, peaks_basename, out_filename, args,
            ratio=args.ratio, only_filter_by='both')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Filter peaks or just add read numbers. \
        Run with either --analysis (-a) and/or --load_bed_info (-b).''')
    #    parser.add_argument('-l', '--load', dest='load',
    #                        default='with_peak_nums/combined_fbf_4_rna_seq.txt',
    #                        help='Load an existing peaks file.')
    parser.add_argument(
        '-d', '--dont_overwrite',
        action='store_true',
        default=False,
        help="""If a given column label already exists, don't overwrite it.""")
    parser.add_argument('-a', '--analysis',
                        default=False,
                        action='store_true',
                        help="Analysis mode: Input file is to be analyzed and filtered (don't add info).")
    parser.add_argument('-i', '--input', dest='input',
                        default='with_peak_nums/',
                        help='Directory of peaks (Required).')
    parser.add_argument('-b', '--load_bed_info',
                        default=False,
                        action='store_true',
                        help="Load bed mode: Add info from bed files to peaks file.")
    parser.add_argument('-o', '--output_dir',
                        default="filtered_<input_dir>/",
                        help='''Output directory for peak files with the reads-in-peak \
                        metrics. This can be the input directory and files will be \
                        edited in place.''')
    parser.add_argument('-r', '--ratio',
                        default='10.',
                        help='''Ratio to filter peaks by.''')
    parser.add_argument('-c', '--config')
    parser.add_argument('-u', '--use',
                        default='fbf1',
                        help='''Which FBF to use (fbf1, fbf2 or both) (Default=fbf1).''')
    parser.add_argument('-f', '--filter_by_both_fbfs',
                         default=False, action='store_true',
                         help='''When use=fbf*_n2 or use=filename,
                         -f will filter both fbfs to the given control.''')
    args = parser.parse_args()
    args.ratio = float(args.ratio)
    import config
    lib = config.config(args.config)
    if args.output_dir == 'filtered_<input_dir>/':
        args.output_dir = 'filtered_' + os.path.dirname(args.input)
    if not os.path.exists(args.output_dir):
        os.system('mkdir ' + args.output_dir)
    if not os.path.exists('tmp_peaks_with_read_nums'):
        os.system('mkdir tmp_peaks_with_read_nums')
    if re.match('.*\.txt', args.input) is not None:
        files_to_process = [args.input]
    else: files_to_process = glob.glob(args.input + '/*.txt')
    if args.load_bed_info:
        start_time = time.time()
        print "Loading .bedgraph info..."
        ga_five = load_all_bedgraph_info(args, lib)
        print "Took %f m to load .bed info." % float((time.time() - start_time)/60.)
    for peaks_filename in files_to_process:
        print "Processing %s" % peaks_filename
        out_filename = args.output_dir + '/' + os.path.basename(peaks_filename)
        peaks = pandas.DataFrame.from_csv(peaks_filename, sep='\t', index_col=False)
        if args.load_bed_info:
            for label in ga_five.keys():
                peaks = add_reads_in_peak(
                    peaks, ga_five[label],
                    dont_overwrite=args.dont_overwrite, label=label)
        peaks.to_csv(out_filename, sep='\t', index_col=False)
        if args.analysis:
            if args.use == 'filename':
                peaks = filter_using_filename(
                    peaks, peaks_filename, out_filename, args)
            elif args.use == 'respective':
                peaks = filter_using_filename(
                    peaks, peaks_filename, out_filename, args)
            else:
                peaks = filter_using_argument(
                    peaks, peaks_filename, out_filename, args)

