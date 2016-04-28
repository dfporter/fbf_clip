import sys
sys.path.insert(0, '/Network/Servers/file.biochem.wisc.edu/Volumes/BioSAN/Users/dfporter/.local/lib/python2.6/site-packages')
sys.path.insert(0, './src/')
#import filter_peaks_by_ratio_in_peak_region
import HTSeq
import pandas
import os
import argparse
import subprocess
#import callpeaks
import matplotlib.pyplot as plt
import matplotlib
from numpy.random import rand
from numpy import arange
import numpy as np
import glob
import re
import shutil
import inspect
import get_sequences
from add_reads_columns import load_bedgraph, add_heights_to_peak_file
#from callpeaks import get_sequences, score_binding_site
#from effect_of_different_controls_figure import *
import add_info_columns
import compare_with_ripchip
import subset_peaks_with_fbe

sys.path.insert(
    0, '/groups/Kimble/Aman Prasad/redo_fbf/analysis/src/feature_locations/'
)
import config

def check_if_all_files_have_read_info(indir):
    for fname in glob.glob(indir + '/*.txt'):
        try:
            peaks = pandas.read_csv(fname, sep='\t')
        except:
            print "Could not open/no rows in %s" % fname
            continue
        cols = set(peaks.columns)
        missing_cols = set(['fbf1_reads', 'fbf2_reads', 'n2_from_fbf1_reads',
                            'n2_from_fbf2_reads', 'rna_seq_reads']) - cols
        if len(list(missing_cols)) > 0:
            print "Missing %s columns in %s" % (
                str(missing_cols), peaks.columns)
            return False
    return True


def add_info_if_needed(indir):
    need_info = False
    for fname in glob.glob(indir + '/*.txt'):
        try:
            peaks = pandas.read_csv(fname, sep='\t')
        except:
            print "Could not open/no reads in %s" % fname
            continue
        cols = set(peaks.columns)
        missing_cols = set(['seq', 'has_fbe']) - cols
        if (len(peaks.index) > 0) and (len(list(missing_cols)) > 0):
            need_info = True
            print "Missing %s columns in %s" % (
                str(missing_cols), peaks.columns)
    if need_info:
        cmd = 'python analysis/src/add_info_columns.py -e -s -d {inputdir}'.format(
            inputdir=indir
        )
        print cmd
        os.system(cmd)


def build_from_raw_peaks():
    cutoff = 5
    for indir in glob.glob('for_methods_comparison/unfiltered/*/'):
        print indir
        #if (indir != 'for_methods_comparison/unfiltered/8_five_reps/') and (
        #    indir != 'for_methods_comparison/unfiltered/9_five_reps/'
        #):
        #    continue
        #do_all_files_have_read_info = check_if_all_files_have_read_info(indir)
        a = '''
    if not do_all_files_have_read_info:
        print "Not all files had read info?"
        ga = {}
        ga['fbf1_reads'] = load_bedgraph('bedgraph_unnorm/combined_fbf1.wig')
        ga['fbf2_reads'] = load_bedgraph('bedgraph_unnorm/combined_fbf2.wig')
        ga['n2_from_fbf1_reads'] = load_bedgraph('bedgraph_unnorm/combined_fbf1_n2.wig')
        ga['n2_from_fbf2_reads'] = load_bedgraph('bedgraph_unnorm/combined_fbf2_n2.wig')
        ga['rna_seq_reads'] = load_bedgraph(
            'bedgraph_unnorm/rna_seq_mod4594.wig')
        cmd = 'python src/filter_peaks_by_ratio_in_peak_region.py -d --load_bed_info '
        cmd += ' -i %s --output_dir %s' % (indir, indir)
        print cmd
        os.system(cmd)'''
    for indir in glob.glob('for_methods_comparison/unfiltered/*/'):
        print "\n*-*-*\nGoing to filter %s" % indir
        #for fname in glob.glob(indir + '/*.txt'):
        #    if not do_all_files_have_read_info:
        #        add_heights_to_peak_file(fname, ga)
        add_info_if_needed(indir)
        for fname in glob.glob(indir + '/*.txt'):
            try:
                peaks = pandas.read_csv(fname, sep='\t')
            except:
                print "Could not open %s" % fname
                continue
           # if len(peaks.index) > 0:
                #if 'biotype' in peaks.columns:
                #    peaks = peaks[peaks['biotype']=='protein_coding']
                #print "No biotype column..."
                #peaks = peaks[peaks['ratio']>=2]
                #peaks.to_csv(fname, sep='\t', index=False)
            peaks = peaks[peaks['ratio']>=cutoff]
            filtered_dir = 'for_methods_comparison/filtered/%s/' % (
                os.path.basename(os.path.dirname(indir)))
            if not os.path.exists('for_methods_comparison/filtered/'):
                os.system('mkdir for_methods_comparison/filtered/')
            if not os.path.exists(filtered_dir):
                os.system('mkdir ' + filtered_dir)
            peaks.to_csv(filtered_dir + '/' + os.path.basename(fname),
                         sep='\t',
                         index=False)
#        filtered_dir = 'for_methods_comparison/filtered/%s/' % os.path.basename(
#            os.path.dirname(indir))
#        cmd = 'python src/filter_peaks_by_ratio_in_peak_region.py --analysis '
#        cmd += '-i {indir_unfilt} --output_dir {outdir_filt} --use both'.format(
#            indir_unfilt=indir, outdir_filt=filtered_dir
#        )
        #os.system(cmd)
    for subdir, dirs, files in os.walk('for_methods_comparison/separate_reps/'):
        for file_basename in files:
            path = os.path.join(subdir, file_basename)
            p = peaks(file=path, name=file_basename)
            fp = p.get_filtered_obj(col='ratio', cutoff=cutoff)
            outname = 'for_methods_comparaison/separate_reps_filtered/'
            dirname = os.path.basename(os.path.dirname(path))
            outname += dirname + '/' + file_basename
            fp.write_table(outname)
    for filtdir in glob.glob('for_methods_comparison/*filtered/*/'):
        print '++\n' * 2
        print filtdir
        add_info_if_needed(filtdir)

'''
Folder fbf1/peaks/combined_fbf1/ contains 2-9 peaks lists.
Folder fbf1/peaks/<rep_name>/null_hyp_*.txt peak lists.
To get to the 5/6 list:
    Copy all 6 fbf*/peaks/<rep_name>/null_hyp_*.txt peak lists to a folder
        (say for_methods_comparison/separate_reps/<number of hypothesis/<rep_name>.txt).
    Run subset...py on this folder and output to another folder
        (say for_methods_comparison/unfiltered/<number_of_hypothesis>/combined_fbf.txt)
    Then filter to:
        (say for_methods_comparison/filtered/<number_of_hypothesis>/combined_fbf.txt)
Add read numbers to these peaks:
    (do we run add_info?)
    python src/filter_peaks*
Run the analysis on these peaks
'''

def stats_for_fig(peaks):
    assert type(peaks) is type([])
    if len(peaks) == 0: return (0,0,0)
    assert type(peaks[0]) is type({})
    known_pos = set(['gld-1', 'htp-1', 'htp-2', 'mpk-1', 'him-3',
                         'fbf-1', 'lip-1', 'syp-2', 'fbf-2', 'fog-1',
                         'fem-3', 'syp-3', 'gld-3', 'fog-3', 'egl-4'])
    num_peaks = len(peaks)
    num_with_fbe = 0.
    found_genes = set()
    for peak in peaks:
        if (peak['has_fbe'] == 1) or (peak['has_fbe'] == '1'):
            num_with_fbe += 1.
        if peak['gene_name'] in known_pos:
            found_genes.add(peak['gene_name'])
    print found_genes
    num_positives = len(list(found_genes))
    perc_fbe = float(100. * num_with_fbe/float(num_peaks))
    return num_peaks, perc_fbe,\
           100. * float(num_positives)/float(len(list(known_pos)))


def make_combined_fig(output_dir='figs/'):
    if not os.path.exists(output_dir):
        os.system('mkdir ' + output_dir)
    # Bar 1: control method 2, no ratio-in-peak control.
    # Bar 2: control method 2, 10 ratio-in_peak control.
    results = {'unfiltered': {}, 'filtered': {}}
    y_tick_labels = []
    sequences = subset_peaks_with_fbe.get_sequences()
    import effect_of_increasing_filter_cutoff as eff
    for control_num in range(2, 10):
        if (control_num == 6) or (control_num == 7): continue
        unfiltered_filename = 'for_methods_comparison/unfiltered/{control_num}_five_reps/combined_fbf.txt'.format(
            control_num=control_num)
        filtered_filename = 'for_methods_comparison/filtered/{control_num}_five_reps/combined_fbf.txt'.format(
            control_num=control_num)
        unfiltered_peaks = pandas.read_csv(unfiltered_filename, sep='\t')
        unfiltered_peaks = subset_peaks_with_fbe.add_seqs(unfiltered_peaks, sequences)
        unfiltered_peaks = subset_peaks_with_fbe.score_binding_site(unfiltered_peaks)
        filtered_peaks = pandas.read_csv(filtered_filename, sep='\t')
        filtered_peaks = subset_peaks_with_fbe.add_seqs(filtered_peaks, sequences)
        filtered_peaks = subset_peaks_with_fbe.score_binding_site(filtered_peaks)
        try:
            print unfiltered_peaks.head(1)
            print filtered_peaks.head(1)
        except:
            print "Error reading %s \n or %s ?" % (unfiltered_filename, filtered_filename)
            continue
        unfiltered_peaks = eff.convert_df_to_list_of_dicts(unfiltered_peaks)
        filtered_peaks = eff.convert_df_to_list_of_dicts(filtered_peaks)
        input_control = 'control_method_' + str(control_num)
        for desc, peaks_df in [('unfiltered', unfiltered_peaks),
                               ('filtered', filtered_peaks)]:
            (num_peaks, perc_fbe, perc_known_pos) = stats_for_fig(peaks_df)
            print perc_known_pos
            y_lab = ''
            if control_num == 2:
                y_lab = 'Local CLIP (Poisson)'
            if control_num == 3:
                y_lab = 'CLIP for gene (Poisson)'
            if control_num == 4:
                y_lab = 'Local control IP (Gauss)'
            if control_num == 5:
                y_lab = 'Control IP for gene (Gauss)'
            if control_num == 6:
                y_lab = 'Local control IP (N.B.)'
            if control_num == 7:
                y_lab = 'Control IP for gene (N.B.)'
            if control_num == 8:
                y_lab = 'Local RNA-seq (Gauss)'
            if control_num == 9:
                y_lab = 'RNA-seq for gene (Gauss)'
            y_lab += ', %s' % desc
            y_tick_labels.append(y_lab)
            results[desc][y_lab] = (num_peaks, perc_fbe, perc_known_pos)
    y_lab_order = [
        'Local CLIP (Poisson), unfiltered',  # 2
        'Local CLIP (Poisson), filtered',
        'CLIP for gene (Poisson), unfiltered',  # 3
        'CLIP for gene (Poisson), filtered',
        'Local control IP (Gauss), unfiltered',  # 4
        'Local control IP (Gauss), filtered',
        'Control IP for gene (Gauss), unfiltered',  # 5
        'Control IP for gene (Gauss), filtered',
        #'Local control IP (N.B.), unfiltered',  # 6
        #'Local control IP (N.B.), filtered',
        #'Control IP for gene (N.B.), unfiltered',  # 7
        #'Control IP for gene (N.B.), filtered',
        'Local RNA-seq (Gauss), unfiltered',  # 8
        'Local RNA-seq (Gauss), filtered',
        'RNA-seq for gene (Gauss), unfiltered',  # 9
        'RNA-seq for gene (Gauss), filtered',
    ]
    print "make_combined_fig(): results['filtered']="
    print results['filtered']
    print "make_combined_fig(): results['unfiltered']="
    print results['unfiltered']
    plot_barchart(
        results, y_lab_order, output_dir=output_dir,
        neg_control='')


def plot_barchart(results, y_labels, output_dir='unknown',
                  neg_control='unknown'):
    '''
            results[desc].append({ylab: (num_peaks, perc_fbe, perc_known_pos)})
    '''
    matplotlib.rc('font', family='Arial')
    if not os.path.exists(output_dir):
        os.system('mkdir ' + output_dir)
    plt.clf()
    _width = 0.35
    ind = np.arange(len(results['filtered']))
    y_labels_filtered = []
    y_labels_unfiltered = []
    y_tick_labels = []
    for n, label in enumerate(y_labels):
        if n % 2:
            y_tick_labels.append(label.partition(',')[0])
        print label
        if re.search(', filtered', label) is not None:
            y_labels_filtered.append(label)
        if re.search(', unfiltered', label) is not None:
            y_labels_unfiltered.append(label)
    pos_1 = ind
    pos_2 = ind + _width
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
    # Plot 1.
    ax1.set_ylabel('Background', fontsize=8)
    unfilt = [results['unfiltered'][x][0] for x in y_labels_unfiltered]
    print unfilt
    ax1.barh(pos_1, unfilt, _width, color='k')
    ax1.barh(pos_2, [results['filtered'][x][0] for x in y_labels_filtered], _width, color='r')
    ax1.set_yticks(ind + _width)
    ax1.set_yticklabels(y_tick_labels)
    ax1.xaxis.grid(True)
    ax1.yaxis.grid(False)
#    ax1.tick_params(axis='y', which='major', bottom='off', top='off')
    max_peak_num = int(sorted([x[0] for x in results['unfiltered'].values()])[-1])
    print 'Ticks put at %s' % str(range(0., max_peak_num + 1000, 1000))
    if max_peak_num > 6000.:
        ax1.set_xticks(
            range(0, max_peak_num + 2000, 2000))
    else:
        ax1.set_xticks(
            range(0, max_peak_num + 1000, 1000))
    ax1.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off')
    ax1.tick_params(axis='y', labelsize=8)
    ax1.axvline(0, color='k', lw=3)
    ax1.set_xlabel('Number of peaks (thousands)', fontsize=8)
    ax1.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: ('%.0f')%(x*1e-3)))
    ax1.set_axisbelow(True)
    no_filt_patch = matplotlib.patches.Patch(color='black', label='None')
    n2_patch = matplotlib.patches.Patch(color='red', label='N2')
    plt.legend(handles=[no_filt_patch, n2_patch],
               markerscale=0.1,
               bbox_to_anchor=(-3.6, .02),
               fontsize=8,
               title='Secondary filter'
               )
    # Plot 2.
    ax2.barh(pos_1, [results['unfiltered'][x][1] for x in y_labels_unfiltered], _width, color='k')
    ax2.barh(pos_2, [results['filtered'][x][1] for x in y_labels_filtered], _width, color='r')
    ax2.axvline(0, color='k', lw=3)
    ax2.set_xlabel('Percent of peaks with FBE', fontsize=8)
    ax2.tick_params(axis='x', labelsize=8)
    ax2.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off')
    ax2.xaxis.grid(True)
    ax2.set_xticks([0, 20, 40, 60, 100], [0, 20, 40, 60, 100])
    ax2.set_axisbelow(True)
    # Plot 3.
    ax3.barh(pos_1, [results['unfiltered'][x][2] for x in y_labels_unfiltered], _width, color='k')
    ax3.barh(pos_2, [results['filtered'][x][2] for x in y_labels_filtered], _width, color='r')
    ax3.set_xticks(range(0, 120, 20))
    ax3.axvline(0, color='k', lw=3)
    ax3.set_xlabel('Percent of known positives identified', fontsize=8)
    ax3.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off')
    ax3.tick_params(axis='x', labelsize=8)
    ax3.xaxis.grid(True)
    ax3.set_axisbelow(True)
    plt.tight_layout(pad=3)#, h_pad=5)
    #plt.subplots_adjust(top=1.5)
    #plt.figtext(0.3, .8, 'Region\nLocal genomic region\nmature RNA')
    #plt.figtext(
    #    0.4, .8, 'Background\nFBF iCLIP\nNo-antibody control iCLIP\nRNA-seq')
    plt.savefig('{subdir}/effect_of_controls_barplot.pdf'.format(
        subdir=output_dir), format='pdf')
    plt.clf()


def complement(s):
    #basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s


def normalize_to_rna_seq(peaks, only_use=['fbf1']):
    fbf1_size =    23367131.
    fbf2_size =    17089670.
    rna_seq_size = 22101373.
    if 'fbf1' in only_use:
        peaks['fbf1_norm_to_rna_seq'] = peaks['fbf1_reads'] * rna_seq_size/fbf1_size
    if 'fbf2' in only_use:
        peaks['fbf2_norm_to_rna_seq'] = peaks['fbf2_reads'] * rna_seq_size/fbf2_size
    return(rna_seq_size/fbf1_size, rna_seq_size/fbf2_size)


def just_info_cols(args):
    rip_targets = compare_with_ripchip.get_ripchip_targets()
    gtf_sep_cols = pandas.read_csv('lib/gtf_with_names_column.txt', sep='\t')
    for control_num in range(2, 8):
        filename = 'with_peak_nums_pypeaks_{control_num}'.format(
            control_num=control_num)
        filename += 'combined_fbf1.txt'
        peaks_df = pandas.read_csv(filename, sep='\t')
        peaks_df = add_info_columns.add_info(peaks_df, rip_targets, gtf_sep_cols)
        peaks_df.to_csv(filename, sep='\t')


def mkdir(_dirname):
    if os.path.isfile(_dirname):
        _dirname = os.path.dirname(_dirname)
    if not os.path.exists(_dirname):
        os.system('mkdir ' + _dirname)


def extract_five_reps_files_from_callpeaks_output_folders(args, lib):
    print inspect.currentframe().f_code.co_name
    rep_folders = []
    for rep_folder in [
        'fbf1_rep_1_all_hyp',
        'fbf1_rep_2_all_hyp',
        'fbf1_rep_3_all_hyp',
        'fbf2_rep_1_all_hyp',
        'fbf2_rep_2_all_hyp',
        'fbf2_rep_3_all_hyp']:
        rep_folders.append(lib[rep_folder])
    print rep_folders
    os.system('rm -r for_methods_comparison/separate_reps/*')
    for null_hyp in range(2, 10):
        fname = 'null_hyp_{i}.txt'.format(i=null_hyp)
        peak_fnames = [rep_folder + '/' + fname for rep_folder in rep_folders]
        print peak_fnames
        reps_under_hyp_dir = lib['for_methods_comparison_separate_reps']
        reps_under_hyp_dir += '/{i}/'.format(i=null_hyp)
        print reps_under_hyp_dir
        if not os.path.exists(reps_under_hyp_dir):
            os.makedirs(os.path.dirname(reps_under_hyp_dir))
        for peak_fname in peak_fnames:
            print '~=~=~'
            print peak_fname
            if not os.path.exists(peak_fname):
                print "%s does not exist." % peak_fname
                continue
            basename = os.path.basename(os.path.dirname(peak_fname)) + '.txt'
            basename = re.sub('run813_fbf', 'fbf2', basename)
            basename = re.sub('runx_fbf', 'fbf1', basename)
            outfile = reps_under_hyp_dir + basename
            print 'copying ' + peak_fname  + ' to ' + outfile
            shutil.copy(peak_fname, outfile)
            peaks = pandas.read_csv(outfile, sep='\t')
            columns =set([x for x in peaks.columns if not (
                    re.search('local', x) or re.search('gene', x))])
            columns = (set(['gene_name']) & set(peaks.columns)) | columns
            columns = columns - set(['exons', 'pvalues'])
            peaks = peaks[list(columns)]
            peaks.to_csv(outfile, sep='\t', columns=list(columns), 
                         index=False)
        cmd = 'python analysis/src/subset_targs_by_all_replicates.py '
        cmd += 'for_methods_comparison/separate_reps/%i/' % null_hyp
        os.system(cmd)
#    for five_reps_folder in glob.glob('for_methods_comparison/separate_reps/*_five_reps/'):
#        shutil.move(five_reps_folder, 'for_methods_comparison/unfiltered/')
    os.system('rm -r {b}/*'.format(
        b='for_methods_comparison/unfiltered/'))
    for five_reps_folder in glob.glob('for_methods_comparison/separate_reps/*_five_reps/'):
        os.system('mv {a} {b}'.format(
            a=five_reps_folder, b='for_methods_comparison/unfiltered/'))
    for filtdir in glob.glob('for_methods_comparison/unfiltered/*/'):
        print '++\n' * 2
        print filtdir
        add_info_if_needed(filtdir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Make a figure of the effect of different controls
on FBE enrichment and peak numbers.
Adds info columns of read numbers in peaks using:
filter_peaks_by_ratio_in_peak_region.py --load_bed_info.

""")
    parser.add_argument('-b', '--build_from_raw_peaks',
                        action='store_true',
                        default=False,
help='Use the raw peaks files (in, say, fbf1/peaks/*) \
 to, for each method of control (2-7) build a 5/6 reps folder, \
 add reads-in-peak data, and create figures of stats vs. control ratio.')
    parser.add_argument('-c', '--combined_fig',
                        action='store_true',
                        default=False,
                        help='Create a combined figure, looking at various \
methods of control. Expects the files output by the --build_from_raw_peaks \
already exist.')
    parser.add_argument('-f', '--start_with_filter',
                        default=False, action='store_true',
                        help='''Use for_methods_comparison/unfiltered/ \
as the input and start by creating for_methods_comparison/filtered/''')
    parser.add_argument('-n', '--info',
                            action='store_true',
                            default=False,
                            help='Just add info columns using \
add_info_columns.py')
    parser.add_argument('-g', '--config', default='config.ini')
    args = parser.parse_args()
    lib = config.config(args.config)
    if args.info:
        just_info_cols(args)
        sys.exit()
    # Get all of the peaks files (all replicates, and the combined)
    # for a given method of control, and put in a folder.
    if args.build_from_raw_peaks:
        print "\t...Building from raw peaks..."
        extract_five_reps_files_from_callpeaks_output_folders(args, lib)
        build_from_raw_peaks()
    if args.start_with_filter:
        build_from_raw_peaks()
    if args.combined_fig:
        mkdir('figs')
        mkdir('methods_comparison')
        # Make the cutoff vs result graph.
        cmd = 'python analysis/src/effect_of_increasing_filter_cutoff.py'
        cmd += ' -i for_methods_comparison/unfiltered/4_five_reps/combined_fbf.txt'
        print cmd
        os.system(cmd)
        make_combined_fig(output_dir='figs/methods_comparison/')
