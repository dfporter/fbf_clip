__author__ = 'dp'
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import build_peaks_arr
import logging
import datetime
import os
logger = logging.getLogger(__name__)

def build_heatmap_from_one_array(
        image, output_filename='heatmap_one_array.pdf', set_scale_from_data=True,
        center_line=True, false_color_part=None):
    from matplotlib import font_manager
    plt.rcParams.update({'font.size': 8})
    original_filename = output_filename
    if image is None or image.size == 0:
        logger.error("Image empty, %s will not be created..." % output_filename)
        return
    x_size = len(image[0])
    narrow_x_width = 100
    print 'false color part: {i}'.format(i=false_color_part)
    if false_color_part is not None:
        for row_index, row in enumerate(image):
            for col_index, val in enumerate(row):
                if false_color_part[row_index][col_index] > 0:
                    image[row_index][col_index] = -0.1 #* image[row_index][col_index]
    for color in ['Reds', 'Spectral', 'coolwarm']:
        for _xrange in ['wide', 'narrow']:
            plt.clf()
            fig, ax1 = plt.subplots()
            plt.xlim([0, x_size])
            plt.ylim([0, len(image)])
            palette = plt.cm.get_cmap(color)
            if set_scale_from_data:
                vmin = np.nanmin(image)
                vmax = np.nanmax(image)
            cax = ax1.imshow(
                image, cmap=palette,
                #extent=[0,len(image[0]), 0, len(image)],
                #aspect=ratio,
                norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=False),
                interpolation='nearest', alpha=0.8)
            if _xrange == 'narrow':
                ax1.set_xlim(x_size/2 - narrow_x_width, x_size/2 + narrow_x_width)
                ax1.set_xticks([x_size/2 - (narrow_x_width),
                                x_size/2 - (narrow_x_width)/2,
                                x_size/2,
                                x_size/2 + narrow_x_width/2,
                                x_size/2 + narrow_x_width])
                ax1.set_xticklabels(['-' + str(narrow_x_width),
                                     '-' + str(narrow_x_width/2),
                                     'Peak center',
                                     str(narrow_x_width/2),
                                     str(narrow_x_width/2)])
            else:
                ax1.set_ylim(0, len(image))
                ax1.set_xticks([0., 500., 1000., 1500., 2000.])
                ax1.set_xticklabels(['-1000', '-500', 'Peak center', '500', '1000'])
            cbar = fig.colorbar(cax, ticks=[vmin, vmax], orientation='horizontal')
            xlab = ['Low', 'High']
            cbar.ax.set_xticklabels(xlab)
            cbar.set_label('Peak height')
            ax1.set_xlabel('Position (nt)', labelpad=20)
            ax1.set_ylabel('Transcript')
            if center_line:
                plt.plot([x_size/2,x_size/2], [0,len(image)], 'k--')
            ax1.set_aspect(1./ax1.get_data_ratio())
            subdir = os.path.dirname(original_filename) + '/{s}/'.format(s=color)
            if not os.path.exists(subdir):
                print 'mkdir ' + subdir
                os.system('mkdir ' + subdir)
            output_filename = subdir + os.path.basename(original_filename).partition('.pdf')[0]
            output_filename += '_%s.pdf' % _xrange
            fig.set_size_inches(4, 4)
            print 'writing to {d}'.format(d=output_filename)
            try:
                plt.savefig(output_filename, format='pdf')
            except:
                logger.error("\t\tCould not save figure %s" % output_filename)
                output_filename = output_filename.partition('.pdf')[0] + '_2.pdf'
                plt.savefig(output_filename, format='pdf')
            logger.info('Output figure %s.' % output_filename)
            plt.clf()
            plt.close()


def build_split_heatmap_from_one_array(
        image_left, image_right, output_filename='figs/region_around_peak.pdf',
        set_scale_from_data=True):
    plt.clf()
    plt.figure(figsize=(6,6))
    plt.axes().set_aspect('equal')
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    palette = plt.cm.get_cmap('Reds')
    palette.set_over('r', 1.0)
    #if set_scale_from_data:
    vmin_l = np.nanmin(image_left)
    vmax_l = np.nanmax(image_left)
    vmin_r = np.nanmin(image_right)
    vmax_r = np.nanmax(image_right)
    ax1.imshow(
        image_left, cmap=palette,
        interpolation='nearest', aspect='auto',
        norm=matplotlib.colors.Normalize(vmin=vmin_l, vmax=vmax_l, clip=False))
    ax1.set_xticks([0., 200., 400., 600., 800., 1000.])
    ax1.set_xticklabels(['-1000', '-800', '-600', '-400', '-200', ''])
    ax1.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='y', labelsize=8)
    ax1.set_xlabel('Distance upstream (nt)', fontdict={'fontsize': 12})
    cax_im = ax2.imshow(
        image_right, cmap=palette,
        interpolation=None,
        #interpolation='nearest',
        aspect='auto',
        norm=matplotlib.colors.Normalize(vmin=vmin_r, vmax=vmax_r, clip=False))
    ax2.set_xticks([0., 200., 400., 600., 800., 1000.])
    ax2.set_xticklabels(['', '200', '400', '600', '800', '1000'])
    ax2.tick_params(axis='x', labelsize=8)
    ax2.tick_params(axis='y', labelsize=8)
    ax2.set_xlabel('Distance downstream (nt)', fontdict={'fontsize': 12})
    #plt.figtext(0.5, 0.05, 'Position relative to peak borders')
    plt.figtext(0.04, 0.5, 'Transcript', rotation=90)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(
        cax_im, cax=cbar_ax,
        ticks=[min(vmin_l, vmin_r), max(vmax_l, vmax_r)],
        orientation='vertical')
    ylab = ['Low', 'High']
    cbar.ax.set_yticklabels(ylab)
    cbar.set_label('CLIP-seq coverage (reads/nt)', fontdict={'fontsize':12})
    try:
        plt.savefig(output_filename, format='pdf')
    except:
        logger.error("\t\tCould not save figure %s" % output_filename)
        output_filename = output_filename.partition('.pdf')[0] + '_2.pdf'
        plt.savefig(output_filename, format='pdf')
    logger.info('Output figure %s.' % output_filename)
    plt.clf()
    plt.close()


def smooth_line(x, y):
    from scipy.interpolate import spline
    import scipy.interpolate
#    xnew = np.linspace(x.min(), x.max(), 300)

#    y_smooth = spline(x, y, xnew)
#    print xnew
#    print y_smooth
    f = scipy.interpolate.interp1d(x, y)#, kind='cubic')
    xnew = np.linspace(x.min(), x.max(), int(1e5))
    plt.clf()
    plt.plot(x, y, 'o', xnew, f(xnew), '--')
    plt.savefig('figs/test_split.pdf', format='pdf')
    plt.clf()
    plt.close()
    return xnew, f(xnew)


def build_heatmap_from_arrays(
        image, image_peaks, output_filename='figs/heatmap.pdf',
        set_scale_from_data=True, include_edge_line=True):
    print "output_heatmap_figures()"
    original_filename = output_filename
    border_pos = []
    if len(image_peaks) < 1:
        print 'print no rows in image.'
        return None
    print image_peaks[0]
    first = isinstance(image_peaks[0][0], np.ma.core.MaskedConstant)
    last = isinstance(image_peaks[0][-1], np.ma.core.MaskedConstant)
    if first and not last:
        aligned_on = 'right'
    elif (not first) and last:
        aligned_on = 'left'
    elif (not first) and (not last):
        first_n = np.isnan(image_peaks[0][0])
        last_n = np.isnan(image_peaks[0][-1])
        if first_n and not last_n:
            aligned_on = 'right'
        elif (not first_n) and last_n:
            aligned_on = 'left'
        else:
            aligned_on = 'middle'
    elif first and last:
        aligned_on = 'middle'
    print aligned_on
    edge_x = []
    if aligned_on == 'left':
        for n, row in enumerate(image_peaks):
            edge_x.append(len(
                [x for x in row if not isinstance(x, np.ma.core.MaskedConstant)])
            )
    if aligned_on == 'right':
        for n, row in enumerate(image_peaks):
            edge_x.append(len(
                [x for x in row if isinstance(x, np.ma.core.MaskedConstant)])
            )
    if aligned_on != 'middle':
        unmasked_x_border, unmasked_y_border = smooth_line(
            np.array(edge_x), np.array(range(0, len(edge_x)))
        )
    # if np.image_peaks[0][0]
    # for row in image_peaks:
    #     for pos in row:
    #         np.isfinite()
    for n, row in enumerate(image_peaks):
        image_peaks[n] = np.ma.masked_where(row <= 0.001, row, copy=True)
    for color in ['Reds', 'Spectral', 'coolwarm']:
        output_filename = original_filename.partition('.pdf')[0] + '_%s.pdf' % color
        plt.clf()
        plt.figure(figsize=(6, 6))
        fig, ax1 = plt.subplots()
        palette = plt.cm.get_cmap(color)
        palette.set_bad('white')
        vmin = np.nanmin(image_peaks)
        vmax = np.nanmax(image_peaks)
        cax = ax1.imshow(
            image_peaks, cmap=palette, #aspect='auto',
            norm=matplotlib.colors.Normalize(
                vmin=vmin, vmax=vmax, clip=False),
            interpolation='nearest', alpha=0.8)
        if include_edge_line and (aligned_on != 'middle'):
            plt.plot(unmasked_x_border, unmasked_y_border, 'k--')
        plt.xlim([
            0, min([len(image_peaks[0]), 5e3])
        ])
        plt.ylim([
            0, min([len(image_peaks), 5e3])
        ])
        plt.gca().invert_yaxis()
        #cbar = fig.colorbar(cax, ticks=[vmin, vmax], orientation='horizontal')
        #xlab = ['Low', 'High']
        #cbar.ax.set_xticklabels(xlab)
        #cbar.set_label('Peak height')
        ax1.set_xlabel('Position (nt)', labelpad=20)
        ax1.set_ylabel('Transcript')
        ax1.set_aspect(1./ax1.get_data_ratio())
        if not os.path.exists(os.path.dirname(output_filename)):
            os.system('mkdir ' + os.path.dirname(output_filename))
        try:
            plt.savefig(output_filename, format='pdf')
            print "Wrote figure to %s." % output_filename
        except:
            logger.error("\t\tCould not save figure %s" % output_filename)
            output_filename = output_filename.partition('.pdf')[0] + '_2.pdf'
            plt.savefig(output_filename, format='pdf')
        logger.info('Output figure %s.' % output_filename)
        plt.clf()
        plt.close()


def plot_intrapeak_distances(peaks, txpts, output_dirname='figs/'):
    logger.info('%s: Plot intrapeak distances called.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    locs = []
    for g in txpts:
        locs.append(len(txpts[g].peak_locs))
    logger.info("average number of peak locs {mn}".format(mn=np.mean(locs)))
    intrapeak_distances = build_peaks_arr.build_intrapeak_distances(peaks, txpts)
    image = []
    output_filename = output_dirname+'/intrapeak_distances.pdf'
    print "*" * 10
    print """Statistics for intrapeak distances:
    Number of distances: {nd}
    Mean distance: {md}
    Median distance: {medd}
    """.format(nd=len(intrapeak_distances),
               md=np.mean(intrapeak_distances),
               medd=np.median(intrapeak_distances))
    max_dist = max(intrapeak_distances)
    for dist in sorted(intrapeak_distances, reverse=True):
        row = [1] * dist
        row += [np.nan] * (max_dist - dist)
        image.append(np.array(row))
    image = np.array(image)
    if image is None or image.size == 0:
        logger.error("Image empty, %s will not be created..." % output_filename)
        return
    plt.clf()
    plt.figure(figsize=(6, 6))
    fig, ax1 = plt.subplots()
    vmin = 0.4
    vmax = 0.5
    palette = plt.cm.get_cmap('Blues')
    palette.set_over('k')
    palette.set_under('white')
    cax = ax1.imshow(
        image, cmap=palette,
        #extent=[0,len(image[0]), 0, len(image)],
        norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=False),
        interpolation='nearest')
    ax1.set_xlim(0, min(len(image[0]), 5000))
    ax1.set_ylim(0, len(image))
    ax1.set_xlabel('Length of intrapeak distance (nt)', labelpad=20)
    ax1.set_ylabel('Peak pair')
    ax1.set_aspect(1./ax1.get_data_ratio())
    try:
        plt.savefig(output_filename, format='pdf')
    except:
        logger.error("\t\tCould not save figure %s" % output_filename)
        output_filename = output_filename.partition('.pdf')[0] + '_2.pdf'
        plt.savefig(output_filename, format='pdf')
    logger.info('Output figure %s.' % output_filename)
    plt.clf()
    plt.close()
