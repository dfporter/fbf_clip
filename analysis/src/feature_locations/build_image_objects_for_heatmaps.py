__author__ = 'dp'
import pandas
import re
import operator
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle
import matplotlib
import copy
import build_peaks_arr
import scipy
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import linkage
import output_heatmap_figures
import logging
import datetime

logger = logging.getLogger(__name__)


def heatmap_of_raw_signal(peaks, txpts, output_dirname='raw_signal_heatmap.pdf',
                          include_motif=False):
    '''Make a heatmap of peak reads in the 3'UTR.
    '''
    if include_motif:
        (utr_dict, signal_in_utr_dict, fbes_dict
         ) = build_peaks_arr.build_peak_raw_utr_signal(
             peaks, txpts, include_motif=include_motif)
    else:
        (utr_dict, signal_in_utr_dict,
         ) = build_peaks_arr.build_peak_raw_utr_signal(
             peaks, txpts, include_motif=include_motif)
    # Signal in utrs is a list of lists.
    utrs = convert_to_list_and_sort(utr_dict, max_len=800)
    peaks = convert_to_list_and_sort(signal_in_utr_dict, max_len=800)
    image_utr_peaks = convert_to_image(peaks, normalize=True, add_padding_to='start', to_absolute_scale=False)
    image_utrs = convert_to_image(utrs, normalize=False, add_padding_to='start', to_absolute_scale=False)
    if include_motif:
        fbes = convert_to_list_and_sort(fbes_dict, max_len=800)
        image_fbes = convert_to_image(fbes, normalize=False, add_padding_to='start', to_absolute_scale=False)
        correlate_signal_and_fbe_location(image_utr_peaks, image_fbes, output_dirname=output_dirname)
    logger.info('%s: Correlate signal and FBE location.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    logger.info('{t}: Creating {f}.'.format(
        t=str(datetime.datetime.now().strftime('%Hh%Mm')),
        f=str(output_dirname+'/raw_signal_in_utrs.pdf')
    ))
    output_heatmap_figures.build_heatmap_from_arrays(
        image_utrs, image_utr_peaks,
        output_filename=output_dirname+'/raw_signal_in_utrs.pdf')
    logger.info('%s: Creating %s.' % (
        datetime.datetime.now().strftime('%Hh%Mm'), str(output_dirname+'/fbes_in_utrs.pdf')
    ))
    if include_motif:
        output_heatmap_figures.build_heatmap_from_arrays(
            image_utrs, image_fbes,
            output_filename=output_dirname+'/fbes_in_utrs.pdf')
    output_heatmap_figures.build_heatmap_from_one_array(
        image_fbes,
        output_filename=output_dirname + '/fbes_in_utrs.pdf')


def correlate_signal_and_fbe_location(coverage, fbes, output_dirname='figs/'):
    sig_at_fbe = {}
    # Are the images aligned on the left or the right?
    if fbes[0].mask[0] and not fbes[0].mask[-1]:
        aligned_on = 'left'
    elif fbes[0].mask[-1] and fbes[0].mask[0]:
        aligned_on = 'right'
    else:
        logger.error("Alignment unclear. Cannot make good correlation.")
        aligned_on = 'left'
    dists = []
    signals = []
    for row_num, row in enumerate(fbes):
        row = row.compressed()
        if not (100 < len(row) < 300): continue
        sig_at_fbe[row_num] = {'dist': [], 'signal': []}
        #pos_of_end_of_utr = find_utr_end(row, aligned_on)
        for pos in range(1, len(row)):
            if row[pos] > 1 and row[pos-1] <= 1:
                sig_at_fbe[row_num]['dist'].append(len(row) - pos)
                sig_at_fbe[row_num]['signal'].append(
                    signal_at_pos_in_row(coverage[row_num].compressed(), pos))
                dists.append(len(row) - pos)
                signals.append(signal_at_pos_in_row(coverage[row_num].compressed(), pos))
    cor_with_fbe_pos = scipy.stats.pearsonr(dists, signals)
    spear_with_fbe_pos = scipy.stats.spearmanr(dists, signals)
    logger.info('Correlation with FBE position: {r} Spearman {s}.'.format(
        r=str(cor_with_fbe_pos), s=str(spear_with_fbe_pos)
    ))
    plt.clf()
    plt.scatter(dists, signals, c='k', alpha=0.5)
    plt.savefig(output_dirname + 'fbe_location_vs_signal.pdf', format='pdf')
    plt.close()
    #sys.exit()


def signal_at_pos_in_row(coverage_row, center_pos):
    left = max(center_pos - 30, 0)
    right = min(center_pos + 30, len(coverage_row) - 1)
    total_coverage = coverage_row[left]
    #for pos in range(left, left + 1): #right):
    #    total_coverage += coverage_row[pos]
    return total_coverage


def find_utr_end(row, aligned_on):
    if aligned_on == 'left':
        for pos, val in enumerate(row):
            if np.isnan(val):
                return pos
    elif aligned_on == 'right':
        for pos, val in enumerate(row):
            if not np.isnan(val):
                return pos
    return False


def heatmap_of_peak_region(peaks, txpts, output_dirname='figs/'):
    """Make a heatmap of peak reads around every peak.
    """
    (all_peaks_raw, region_around_all_peaks_raw, highest_peak_raw,
     secondary_peaks_raw, region_around_all_peaks_raw,
     all_peaks_highest_marked, secondary_peaks_highest_marked,
     highest_if_have_secondary_raw,
     highest_if_have_secondary_highest_marked
     ) = build_peaks_arr.build_peak_raw_peak_region(peaks, txpts)
    # Non-split heatmaps.
    # All peaks.
    image = convert_to_image(all_peaks_raw)
    all_peaks_highest_marked = convert_to_image(all_peaks_highest_marked)
    secondary_peaks_highest_marked = convert_to_image(secondary_peaks_highest_marked)
    # Image.
    logger.info('%s: Creating all_peaks_hm_by_height.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    output_heatmap_figures.build_heatmap_from_one_array(
        image,
        false_color_part=all_peaks_highest_marked,
        output_filename=output_dirname + '/all_peaks_hm_by_height.pdf', center_line=True)
    # Image.
    logger.info('%s: Creating highest_peak_hm_by_height.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    image = convert_to_image(highest_peak_raw)
    output_heatmap_figures.build_heatmap_from_one_array(
        image, output_filename=output_dirname + '/highest_peak_hm_by_height.pdf', center_line=True)
    # Image.
    logger.info('%s: Creating secondary_peaks_hm_by_height.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    if len(secondary_peaks_raw) > 0:
        image = convert_to_image(secondary_peaks_raw)
        output_heatmap_figures.build_heatmap_from_one_array(
            image,
            false_color_part=secondary_peaks_highest_marked,
            output_filename=output_dirname + '/secondary_peaks_hm_by_height.pdf', center_line=True)
    #image_left = convert_to_image(region_around_all_peaks_raw[0])
    #image_right = convert_to_image(region_around_all_peaks_raw[1])
    #logger.info('%s: Creating region_around_peak_hm_by_height.' % (
    #    datetime.datetime.now().strftime('%Hh%Mm')))
    #output_heatmap_figures.build_split_heatmap_from_one_array(
    #    image_left, image_right,
    #    output_filename=output_dirname + '/region_around_peak_hm_by_height.pdf',
    #    set_scale_from_data=True)
    heatmaps_by_similarity(
        all_peaks_raw, highest_peak_raw, secondary_peaks_raw,
        region_around_all_peaks_raw,
        all_peaks_highest_marked,
        secondary_peaks_highest_marked,
        highest_if_have_secondary_raw,
        highest_if_have_secondary_highest_marked,
        metric='correlation',
        output_dirname=output_dirname)


def heatmaps_by_similarity(
        all_peaks_raw, highest_peak_raw, secondary_peaks_raw,
        region_around_all_peaks_raw,
        all_peaks_highest_marked,
        secondary_peaks_highest_marked,
        highest_if_have_secondary_raw,
        highest_if_have_secondary_highest_marked,
        metric='correlation',
        output_dirname='figs/'):
    # Image. Highest if secondary in range, sorted by proximity.
    sub_i, sub_m = subset_to_secodary_with_proximal_primary(
        highest_if_have_secondary_raw,
        highest_if_have_secondary_highest_marked
    )
    print ' Highest if secondary in range, sorted by proximity.'
    print type(sub_i)
    print type(sub_m)
    print sub_i
    print sub_m
    (sub_i, sub_m) = cluster_image_by_distance_to_mask(
        convert_to_image(sub_i), second_image=sub_m)
    output_heatmap_figures.build_heatmap_from_one_array(
        sub_i, #false_color_part=highest_if_have_secondary_highest_marked,
        output_filename=output_dirname + '/highest_peak_if_have_secondary_in_range_peak_highest_marked_hm_by_proximity.pdf', center_line=True)
    # Image. Highest if secondary in range, sorted by correlation.
    sub_i, sub_m = subset_to_secodary_with_proximal_primary(
        highest_if_have_secondary_raw,
        highest_if_have_secondary_highest_marked
    )
    tup = cluster_image(
        convert_to_image(sub_i), false_color_part=sub_m,
        metric=metric)
    print tup
    if len(tup) == 2:
        (sub_i, sub_m) = tup
        output_heatmap_figures.build_heatmap_from_one_array(
            sub_i, #false_color_part=highest_if_have_secondary_highest_marked,
            output_filename=output_dirname + '/highest_peak_if_have_secondary_in_range_peak_highest_marked_hm_by_%s.pdf' % metric, center_line=True)
        # Image.
        sub_i, sub_m = subset_to_secodary_with_proximal_primary(
            secondary_peaks_raw,
            secondary_peaks_highest_marked
        )
    tup = cluster_image(
        convert_to_image(sub_i),
        metric=metric, false_color_part=sub_m)
    print tup
    if len(tup) == 2:
        sub_i, sub_m = tup
        output_heatmap_figures.build_heatmap_from_one_array(
            sub_i, false_color_part=sub_m,
            output_filename=output_dirname + '/secondary_peaks_highest_marked_only_those_with_primary_in_range_hm_by_%s.pdf' % metric, center_line=True)
    # Image.
    logger.info('%s: Creating all_peaks_hm_by*.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    logger.info('%s:\tConverting to image.' % (datetime.datetime.now().strftime('%Hh%Mm')))
    image = convert_to_image(all_peaks_raw)
    all_peaks_highest_marked = convert_to_image([list(x) for x in all_peaks_highest_marked])
    secondary_peaks_highest_marked = convert_to_image([list(x) for x in secondary_peaks_highest_marked])
    logger.info('%s:\tClustering.' % (datetime.datetime.now().strftime('%Hh%Mm')))
    (image, all_peaks_highest_marked) = cluster_image(
        image, metric=metric, false_color_part=all_peaks_highest_marked)
    logger.info('%s:\tWriting image.' % (datetime.datetime.now().strftime('%Hh%Mm')))
    output_heatmap_figures.build_heatmap_from_one_array(
        image, false_color_part=all_peaks_highest_marked,
        output_filename=output_dirname + '/all_peaks_hm_by_%s.pdf' % metric, center_line=True)
    # Image.
    logger.info('%s: Creating highest_peak_hm_by_*.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    output_heatmap_figures.build_heatmap_from_one_array(
        cluster_image(convert_to_image(highest_peak_raw), metric=metric),
        output_filename=output_dirname + '/highest_peak_hm_by_%s.pdf' % metric, center_line=True)
    # Image.
    logger.info('%s: Creating highest_peak_if_have_secondary_peak_hm_by_*.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    output_heatmap_figures.build_heatmap_from_one_array(
        cluster_image(convert_to_image(highest_if_have_secondary_raw), metric=metric),
        output_filename=output_dirname + '/highest_peak_if_have_secondary_peak_hm_by_%s.pdf' % metric, center_line=True)
    #print highest_if_have_secondary_highest_marked
    # Image.
    if len(highest_if_have_secondary_raw) > 0 and len(highest_if_have_secondary_highest_marked)>0:
        (image, highest_if_have_secondary_highest_marked) = cluster_image(
            convert_to_image(highest_if_have_secondary_raw), metric=metric,
            false_color_part=highest_if_have_secondary_highest_marked)
        output_heatmap_figures.build_heatmap_from_one_array(
            image, false_color_part=highest_if_have_secondary_highest_marked,
            output_filename=output_dirname + '/highest_peak_if_have_secondary_peak_highest_marked_hm_by_%s.pdf' % metric, center_line=True)
    # Image.
    if len(secondary_peaks_raw) > 0 and len(secondary_peaks_highest_marked) > 0:
        (image, secondary_peaks_highest_marked) = cluster_image(
            convert_to_image(secondary_peaks_raw),
            metric=metric, false_color_part=secondary_peaks_highest_marked)
        output_heatmap_figures.build_heatmap_from_one_array(
            image,
            output_filename=output_dirname + '/secondary_peaks_hm_by_%s.pdf' % metric,
            center_line=True)
        # Image.
        output_heatmap_figures.build_heatmap_from_one_array(
            image, false_color_part=secondary_peaks_highest_marked,
            output_filename=output_dirname + '/secondary_peaks_highest_marked_hm_by_%s.pdf' % metric, center_line=True)
    else: print "Number of secondary peaks {i}. Numbe of secondary peaks with the highest marked {i2}.".format(
        i=len(secondary_peaks_raw),
        i2=len(secondary_peaks_highest_marked)
    )
    # Image. This one is currently broken.
#    image_left = convert_to_image(region_around_all_peaks_raw[0])
#    image_right = convert_to_image(region_around_all_peaks_raw[1])
#    ordered_image1, ordered_image2 = cluster_images(image_left, image_right, metric=metric)
#    output_heatmap_figures.build_split_heatmap_from_one_array(
#        ordered_image1, ordered_image2,
#        output_filename=output_dirname + '/region_around_peak_hm_by_%s.pdf' % metric,
#        set_scale_from_data=True)


def subset_to_secodary_with_proximal_primary(image, mask):
    logger.info('Subsetting image to those rows with a nonzero mask value.')
    logger.info('Image is {r} rows and mask is {m} rows.'.format(
        r=len(image), m=len(mask)
    ))
    new_image = []
    new_mask = []
    for index, row in enumerate(image):
        if np.max(mask[index]) > 0:
            new_image.append(row)
            new_mask.append(mask[index])
    logger.info('Subsetted image is {r} rows.'.format(r=len(new_image)))
    return np.array(new_image), np.array(new_mask)


def cluster_image_by_distance_to_mask(image, second_image=None):
    check_type(image)
    if len(image) == 0:
        if second_image is None:
            return image
        return image, second_image
    check_type(image)
    middle = int(len(image[0])/2.)
    index_to_closest_dist = {}
    for row_n, row in enumerate(image):
        if np.max(second_image[row_n]) <= 0:
            continue
        borders = []
        last_val = 0
        for index, val in enumerate(second_image[row_n]):
            if (val > 0) and (last_val <= 0):
                if len(borders) == 0:
                    borders.append({'left': index})
                    continue
                elif 'right' not in borders[-1]:  # Just replace if the right border missing.
                    borders[-1] = {'left': index}
                else:
                    borders.append({'left': index})
            if (val <= 0) and (last_val > 0):
                borders[-1]['right'] = index
            last_val = val
        if 'right' not in borders[-1]:
            borders[-1]['right'] = len(image[n]-1)
        index_to_closest_dist[row_n] = 1e4
        for peak_range in borders:
            if len(peak_range.keys()) != 2: continue
            dist = np.min([
                abs(peak_range['right']-middle), abs(peak_range['left']-middle)])
            index_to_closest_dist[row_n] = np.min([index_to_closest_dist[row_n], dist])
    print 'presort'
    print index_to_closest_dist
    print 'keys'
    print index_to_closest_dist.keys()
#    index_list = [x for x in index_to_closest_dist.keys()]
    sorted_indexes = sorted(
        index_to_closest_dist.keys(), reverse=True,
        key=lambda x: index_to_closest_dist[x])
    print sorted_indexes
    image = image[sorted_indexes,]
    second_image = second_image[sorted_indexes,]
    return image, second_image

def check_type(image):
    try:
        assert (type(image) in [type(np.array([])), type(np.ma.core.MaskedArray([]))])
    except:
        print "unexpected type: %s" % type(image)

def cluster_image(
        image, metric='correlation', do_second_image=False, second_image=None,
        false_color_part=None):
    check_type(image)
    if len(image) == 0:
        if second_image is None:
            return image
        return image, second_image
    check_type(image[0])
    dists_image = scipy.spatial.distance.pdist(image, metric=metric)
    Z = sch.linkage(dists_image)
    heatmap_order = sch.leaves_list(Z)
    print 'cluster image: image type:{i} false color type {f}'.format(
        i=type(image), f=type(false_color_part)
    )
    if type(false_color_part) == np.ndarray:
        print "Types of the first row: image {i} false color {f}".format(
            i=type(image[0]), f=type(image[0])
        )
        print 'size image {i}, size false color {f}, ' \
              'colums image {c} columns false color {fc}'.format(
            i=len(image), f=len(false_color_part),
            c=len(image[0]), fc=len(false_color_part[0])
        )
    ordered_image = image[heatmap_order,:]
    if false_color_part is not None:
        false_color_part = false_color_part[heatmap_order,:]
    if not do_second_image:
        if false_color_part is not None: return (ordered_image, false_color_part)
        else: return ordered_image
    if do_second_image:
        second_image = second_image[heatmap_order,:]
        if false_color_part is not None:
            return (ordered_image, second_image, false_color_part)
        else:
            return (ordered_image, second_image)


def cluster_images(image1, image2, metric='correlation'):
    check_type(image1)
    check_type(image2)
    new_image = np.concatenate((image1, image2), axis=1)
    dists_image = scipy.spatial.distance.pdist(new_image, metric=metric)
    Z = sch.linkage(dists_image)
    heatmap_order = sch.leaves_list(Z)
    ordered_image1 = image1[heatmap_order,:]
    ordered_image2 = image2[heatmap_order,:]
    return ordered_image1, ordered_image2


def convert_to_image(
        list_of_lists, normalize=True, add_padding_to='end', to_absolute_scale=False):
    if len(list_of_lists) == 0:
        print 'empty cols'
        return np.array(np.array([]))
    for index, alist in enumerate(list_of_lists):
        if len(alist) < 2:
            del list_of_lists[index]
    longest_loc = np.max([len(t) for t in list_of_lists])
    print 'convert_to_image:'
    print longest_loc
    padded_lists = pad_list(list_of_lists, longest_loc, add_to=add_padding_to)
    filtered_lists = []
#    for index, row in enumerate(padded_lists):
#        if len(row[~np.isnan(row)]) > 0:
#            filtered_lists.append(row)
    padded_arrays = np.array(padded_lists)
    print 'non-nan'
    print padded_arrays[0][~np.isnan(padded_arrays[0])]
    print 'w nan'
    print padded_arrays[0]
    arr_2d_image = np.array(padded_arrays)  # image is now a 2d np.array.
    image = padded_arrays
    if normalize:
        image = normalize_rows(padded_arrays, to_absolute_scale=to_absolute_scale)
    print 'image:'
    print image
    image = np.ma.masked_where(np.isnan(image), image)
    print 'after mask'
    print image
    return image


def heatmap_by_gene_length_just_peak_ranges(peaks, txpts, output_dirname='figs/'):
    """Sort targets by length and output a heatmap of binding site locations.
    peaks : pandas.DataFrame
    txpts : dict by gene name of flocs objects
    """
    (locsarr, peakslocs, utrsarr, peaks_in_utrs) = build_peaks_arr.build_peak_dicts(peaks, txpts)
    (padded_list, padded_peaks, padded_utrs, padded_peaks_in_utrs) = sort_by_length_and_pad(
        locsarr, peakslocs, utrsarr, peaks_in_utrs)
    #masked_padded = np.array([np.ma.array(t, mask=np.isnan(t)) for t in padded_peaks])
    image_txpt = np.array([t for t in padded_list])
    image_peaks = np.array([t for t in padded_peaks])
    image_txpt = np.ma.masked_where(np.isnan(image_txpt), image_txpt)
    image_peaks = np.ma.masked_where(np.isnan(image_peaks), image_peaks)
    output_heatmap_figures.build_heatmap_from_arrays(
        image_txpt, image_peaks, output_filename=output_dirname + '/peaks_in_txpt_heatmap.pdf')
    # Create a map of the 3'UTRs.
    image_utrs = np.array([t for t in padded_utrs])
    image_utr_peaks = np.array([t for t in padded_peaks_in_utrs])
    output_heatmap_figures.build_heatmap_from_arrays(
        image_utrs, image_utr_peaks, output_filename=output_dirname + '/peaks_in_utrs_heatmap.pdf')


def normalize_rows(image, to_absolute_scale=False):
    '''Convert a 2d list/array of lists/arrays to an array of arrays.
    '''
    norm_rows = []
    for row in image:
        #_max = np.nanmax([1., float(np.nanmax(row))])
        _max = np.nanmax(row)
        if to_absolute_scale:
            norm_rows.append(np.array(
                [min(t, to_absolute_scale) for t in row]
            ))
        else:
            if np.isnan(_max) or _max <= 1:
                norm_rows.append(row)
            else:
                norm_rows.append(np.array([t/_max for t in row]))
    return np.array(norm_rows)


def pad_list(locs, pad_to, pad_with="NA", add_to='end'):
    """Pad a list of lists/arrays and return a list of lists.
    """
    padded_list = []
    for t in locs:
        if len(t) == 0: continue
        if len(t) < pad_to:
            if add_to == 'end':
                padded_list.append(t + [pad_with] * (pad_to-len(t)))
            elif add_to == 'start':
                padded_list.append([pad_with] * (pad_to-len(t)) + t)
        else:
            padded_list.append(t)
    for row in padded_list:
        if row is None:
            logger.warn("empty row?")
            continue
        for index, item in enumerate(row):
            if item == 'NA':
                row[index] = np.nan
    return padded_list


def convert_to_list_and_sort(arr, max_len=5000):
    '''Converts a list/array of lists/arrays to a list of lists/arrays, sorted.
    '''
    locslist = [arr[x] for x in arr if (10 < len(arr[x]) < max_len)]
    locslist.sort(key=lambda t: len(t))
    return locslist


def sort_by_length_and_pad(locsarr, peakslocs, utrsarr, peaks_in_utrs):
    locslist = convert_to_list_and_sort(locsarr)
    peakslocs = convert_to_list_and_sort(peakslocs)
    utrslist = convert_to_list_and_sort(utrsarr)
    peaks_in_utrs_list = convert_to_list_and_sort(peaks_in_utrs)
#    locslist = [locsarr[x] for x in locsarr if (len(locsarr[x]) < 5000)]
#    peakslocs = [peakslocs[x] for x in peakslocs if len(peakslocs[x]) < 5000]
#    locslist.sort(key=lambda t: len(t))
#    peakslocs.sort(key=lambda t: len(t))
    longest_loc = np.max([len(t) for t in locslist])
    padded_list = pad_list(locslist, longest_loc)
    padded_peaks = pad_list(peakslocs, longest_loc)
    padded_utrs = pad_list(utrslist, np.max([len(t) for t in utrslist]))
    padded_utr_peaks = pad_list(peaks_in_utrs_list, np.max([len(t) for t in utrslist]))
    return padded_list, padded_peaks, padded_utrs, padded_utr_peaks


