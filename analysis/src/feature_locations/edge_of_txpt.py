import build_image_objects_for_heatmaps
import output_heatmap_figures
import numpy as np
import matplotlib.pyplot as plt
import os


def clip_end_of_txpt(list2d, edge=200, which_end='left'):
    end_of_list2d = []
    for row in list2d:
        if len(row) >= (edge + 1):
            if which_end == 'right':
                end_of_list2d.append(row[-edge:])
            if which_end == 'left':
                end_of_list2d.append(row[:edge])
    return end_of_list2d


def make_graph_of_a_txpt_region(list2d, list2d_fbe, edge=200, outfname='figs/200.pdf'):
    print "Making a graph of a transcript region (5' and 3' %i nt)." % int(edge)
    if not os.path.exists(os.path.dirname(outfname)):
        os.system('mkdir ' + os.path.dirname(outfname))
    make_left_graph(list2d, list2d_fbe, edge=edge,
                    outfname=outfname.partition('.pdf')[0] + '_left.pdf')
    make_right_graph(list2d, list2d_fbe, edge=edge,
                    outfname=outfname.partition('.pdf')[0] + '_right.pdf')


def make_left_graph(list2d, list2d_fbe, edge=200, outfname='figs/200.pdf'):
    end_of_list2d = clip_end_of_txpt(list2d, edge=edge,
                                     which_end='left')
    end_of_list2d_fbe = clip_end_of_txpt(list2d_fbe, edge=edge,
                                         which_end='left')
    end_of_image = build_image_objects_for_heatmaps.convert_to_image(end_of_list2d)
    end_of_image_fbe = build_image_objects_for_heatmaps.convert_to_image(end_of_list2d_fbe, normalize=False)
    output_heatmap_figures.build_heatmap_from_arrays(
        None, end_of_image,
        output_filename=outfname.partition('.pdf')[0] + '_200.pdf',
        set_scale_from_data=True)
    output_heatmap_figures.build_heatmap_from_arrays(
        None, end_of_image_fbe,
        output_filename=outfname.partition('.pdf')[0] + '_200_fbe.pdf',
        set_scale_from_data=True)
    plot_a_couple_aves(end_of_image, end_of_image_fbe, outfname,
                       which_end='left')
    end_of_image = build_image_objects_for_heatmaps.cluster_image(
        end_of_image,
        do_second_image=False,
        metric='euclidean', false_color_part=None)
    end_of_image_fbe = build_image_objects_for_heatmaps.cluster_image(
        end_of_image_fbe,
        do_second_image=False,
        metric='euclidean', false_color_part=None)
    output_heatmap_figures.build_heatmap_from_arrays(
        None, end_of_image, include_edge_line=False,
        output_filename=outfname.partition('.pdf')[0] + '_200_clustered.pdf',
        set_scale_from_data=True)
    output_heatmap_figures.build_heatmap_from_arrays(
        None, end_of_image_fbe, include_edge_line=False,
        output_filename=outfname.partition('.pdf')[0] + '_200_clustered_fbe.pdf',
        set_scale_from_data=True)


def make_right_graph(list2d, list2d_fbe, edge=200, outfname='figs/200.pdf'):
    end_of_list2d = clip_end_of_txpt(list2d, edge=edge,
                                     which_end='right')
    end_of_list2d_fbe = clip_end_of_txpt(list2d_fbe, edge=edge,
                                         which_end='right')
    end_of_image = build_image_objects_for_heatmaps.convert_to_image(end_of_list2d)
    end_of_image_fbe = build_image_objects_for_heatmaps.convert_to_image(end_of_list2d_fbe, normalize=False)
    output_heatmap_figures.build_heatmap_from_arrays(
        None, end_of_image, include_edge_line=False,
        output_filename=outfname.partition('.pdf')[0] + '_200.pdf',
        set_scale_from_data=True)
    output_heatmap_figures.build_heatmap_from_arrays(
        None, end_of_image_fbe, include_edge_line=False,
        output_filename=outfname.partition('.pdf')[0] + '_200_fbe.pdf',
        set_scale_from_data=True)
    plot_a_couple_aves(end_of_image, end_of_image_fbe, outfname,
                       which_end='right')
    end_of_image = build_image_objects_for_heatmaps.cluster_image(
        end_of_image,
        do_second_image=False,
        metric='euclidean', false_color_part=None)
    end_of_image_fbe = build_image_objects_for_heatmaps.cluster_image(
        end_of_image_fbe,
        do_second_image=False,
        metric='euclidean', false_color_part=None)
    output_heatmap_figures.build_heatmap_from_arrays(
        None, end_of_image, include_edge_line=False,
        output_filename=outfname.partition('.pdf')[0] + '_200_clustered.pdf',
        set_scale_from_data=True)
    output_heatmap_figures.build_heatmap_from_arrays(
        None, end_of_image_fbe, include_edge_line=False,
        output_filename=outfname.partition('.pdf')[0] + '_200_clustered_fbe.pdf',
        set_scale_from_data=True)


def plot_a_couple_aves(end_of_image, end_of_image_fbe, outfname,
                       which_end='left'):
    aves = np.average(end_of_image, axis=0)
    aves_fbe = np.average(end_of_image_fbe, axis=0)
    ave_first_half = np.average(aves[:100])
    ave_first_half_fbe = np.average(aves_fbe[:100])
    plt.clf()
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(range(len(aves)), aves, 'k-'),
    ax1.set_xlabel('Nucleotides')
    ax1.set_ylabel('Average FBF iCLIP coverage\n(fraction of maximum density)')
    ax2.plot(range(len(aves_fbe)), aves_fbe, 'r-')
    ax2.set_xlabel('Nucleotides')
    ax2.set_ylabel('FBE density\n(average probability of an FBE)')
    if which_end == 'left':
        ax1.set_xticklabels(["5' end", '50', '100', '150', '200'])
        ax2.set_xticklabels(["5' end", '50', '100', '150', '200'])
    if which_end == 'right':
        ax1.set_xticklabels(["-200", '-150', '-100', '-50', "3' end"])
        ax2.set_xticklabels(["-200", '-150', '-100', '-50', "3' end"])
    plt.savefig(outfname.partition('.pdf')[0] + '_averages.pdf')
    plt.clf()
    plt.close()

