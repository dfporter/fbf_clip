"""
Test the following, quoted from the Zhang lab wobsite:
<quote>
Since CLIP tags with crosslinking induced mutations (i.e. deletions)
most likely represent tags that read through the cross link site, it
is recommended that these tags are removed.

cut -f 4 tag.uniq.del.bed | sort | uniq > tag.uniq.del.id
python ~/src/CIMS/joinWrapper.py tag.uniq.trunc.bed tag.uniq.del.id \
    4 1 V tag.uniq.trunc.clean.bed

Here tag.uniq.del.bed is the list of deletions in unique tags.
You should have obtained this file by CIMS analysis. 
<\quote>
"""
import glob
import sys
import os
import re


def remove_deletions(shifted_filename='', do_all=False):
    if not os.path.exists('no_deletions/'):
        os.system('mkdir no_deletions')
    if not os.path.exists('temp_dir/'):
        os.system('mkdir temp_dir')
    if do_all:
        for input_filename in glob.glob('shifted_novo_tags/*.bed'):
            remove_deletions_this_file(input_filename)
    elif shifted_filename != '':
        remove_deletions_this_file(shifted_filename)


def remove_deletions_this_file(input_filename):
    dels_filename = 'mismatches_by_type/' + os.path.basename(
        input_filename).partition('.bed')[0] + '_del.bed'
    if not os.path.exists(dels_filename):
        print "No deletions file; expected: %s" % dels_filename
        return
#    temp_uniq = 'temp_dir/' + os.path.basename(input_filename).partition(
#        '.bed')[0] + '_temp_unique'
    tag_ids_filename = 'temp_dir/' + os.path.basename(input_filename).partition(
        '.bed')[0] + '_tag_ids'
    cmd = 'cut -f 4 {del_tags} | sort | uniq > {tag_ids}'.format(
        del_tags=dels_filename, 
        tag_ids=tag_ids_filename)
    print cmd
    os.system(cmd)
    # joinWrapper command.
    clean_filename = 'no_deletions/' + os.path.basename(input_filename)
    cmd = 'python ../clip/CIMS/joinWrapper.py '
    cmd += '{trunc} {del_id} 4 1 V {clean}'.format(
        trunc=input_filename,
        del_id=tag_ids_filename,
        clean=clean_filename)
    print cmd
    os.system(cmd)


if __name__ == '__main__':
    remove_deletions(do_all=True)



