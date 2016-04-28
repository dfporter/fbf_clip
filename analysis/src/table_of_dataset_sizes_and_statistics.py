import sys
import os
import subprocess
import re
import xlwt
import argparse
import glob
sys.path.insert(
    0, '/groups/Kimble/Aman Prasad/redo_fbf/analysis/src/feature_locations/')
import config


def dict_of_number_of_lines(_list, cofactor=1):
    return dict([(x, cofactor * len(open(x).readlines())) for x in _list])

def fastq_stats(lib):
    if 'fbf1_fastq' in lib:
        fbf1_fastq = lib['fbf1_fastq']
    else:
        fbf1_fastq = 'fastq/fbf1/'
    if 'fbf2_fastq' in lib:
        fbf2_fastq = lib['fbf2_fastq']
    else:
        fbf2_fastq = 'fastq/fbf2/'
    fbf1_exp_files, fbf1_n2_files = load_dir(fbf1_fastq + '/*fastq')
    fbf2_exp_files, fbf2_n2_files = load_dir(fbf2_fastq + '/*fastq')
    li = 'Reads in fastq files:'
    li += write_lines(fbf1_exp_files, fbf1_n2_files, fbf2_exp_files,
                fbf2_n2_files, cofactor=0.25)
    print li
    return li

def write_lines(fbf1_exp_files, fbf1_n2_files, fbf2_exp_files,
                fbf2_n2_files, cofactor=1.):
    li = '''
    FBF-1 exp files:
    {e1}
    FBF-1 N2 files:
    {n1}
    FBF-2 exp files:
    {e2}
    FBF-2 N2 files:
    {n2}
'''.format(
    e1=dict_of_number_of_lines(fbf1_exp_files, cofactor=cofactor),
    n1=dict_of_number_of_lines(fbf1_n2_files, cofactor=cofactor),
    e2=dict_of_number_of_lines(fbf2_exp_files, cofactor=cofactor),
    n2=dict_of_number_of_lines(fbf2_n2_files, cofactor=cofactor))
    return li

def load_dir(_fastq):
    exp_files = [x for x in glob.glob(_fastq) if \
                    re.search('exp', x)]
    n2_files = [x for x in glob.glob(_fastq) if \
                    re.search('n2', x)]
    return exp_files, n2_files

def bed_stats():
    fbf1, _ = load_dir('fbf1/bed_collapsed/*.bed')
    _, fbf1_n2 = load_dir('new_n2/bed_collapsed/*fbf1*.bed')
    fbf2, _ = load_dir('fbf2/bed_collapsed/*.bed')
    _, fbf2_n2 = load_dir('new_n2/bed_collapsed/*fbf2*.bed')
    li = "Reads in collapsed bed files:"
    li += write_lines(fbf1, fbf1_n2, fbf2, fbf2_n2)
    print li
    return li

def write_workbook(stats):
    book = xlwt.Workbook()
    sh = book.add_sheet('stats')
    sh.write(0, 0, 'Dataset')
    sh.write(0, 1, 'Total reads')
    sh.write(0, 2, 'Mapped')
    sh.write(0, 3, '% Mapped')
    sh.write(0, 4, 'Mapped at 20 MAPQ')
    n = 1
    for filename in stats:
        sh.write(n, 0, filename)
        sh.write(n, 1, stats[filename]['total'])
        sh.write(n, 2, stats[filename]['mapped'])
        if stats[filename]['total'] == 0:
            perc_mapped = 0
        else:
            perc_mapped = float(stats[filename]['mapped'])/float(stats[filename]['total'])
        sh.write(n, 3, '%.3f' % float(100.0 * perc_mapped))
        sh.write(n, 4, stats_mapq[filename]['mapped'])
        n += 1
    book.save('./tables/mapping_quality_stats.xls')
    #bam_filename = sys.argv[1]
    #print dataset_size(bam_filename)


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_ini', default='analysis/src/config.ini')
#parser.add_argument(
args = parser.parse_args()
lib = config.config(filepath=args.config_ini)
#fastq_stats(lib)
bed_stats()
stats_mapq = {}
stats = {}
#li = "Dataset\tTotal reads\tMapped\tMapped at 20 MAPQ\n"
#for filename in stats:
#    li += '{name}\t{total}\t{mapped}\t{mapq}\n'.format(
#        name=filename, total=stats[filename]['total'],
#        mapped=stats[filename]['mapped'], mapq=stats_mapq[filename]['mapped'])
#print li
#with open('./tables/mapping_quality_stats.txt', 'w') as f:
#    f.write(li)
