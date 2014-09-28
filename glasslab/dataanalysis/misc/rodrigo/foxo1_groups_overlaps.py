'''
Created on Sep 28, 2014

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.dataanalysis.misc.rodrigo.samples import sample_name,\
    get_threshold
if __name__ == '__main__':
    yzer = SeqGrapher()

    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/' +\
        'Miscellaneous_Collaborations/Rodrigo_CD8s_2014_09/Enhancers'
    dirpath = yzer.get_path(dirpath)

    cond, seq, breed = ('naive', 'atac', '')
    wt_prefix = sample_name(cond, seq, breed)
    ko_prefix = sample_name(cond, seq, 'foxo1_ko_')

    wt_dirpath = yzer.get_filename(dirpath, wt_prefix)
    ko_dirpath = yzer.get_filename(dirpath, ko_prefix)

    wt_filename = yzer.get_filename(wt_dirpath,
                                    wt_prefix + '_enhancers.txt')
    ko_filename = yzer.get_filename(ko_dirpath,
                                    ko_prefix + '_enhancers.txt')

    wt_data = yzer.import_file(wt_filename)
    wt_data = wt_data.fillna(0)
    ko_data = yzer.import_file(ko_filename)
    ko_data = ko_data.fillna(0)

    min_thresh = get_threshold(seq)
    wt_data = wt_data[wt_data['tag_count'] >= min_thresh]
    ko_data = ko_data[ko_data['tag_count'] >= min_thresh]

    wt_only = wt_data[
        wt_data['foxo1_ko_naive_atac_tag_count'] < min_thresh]

    both = wt_data[
        wt_data['foxo1_ko_naive_atac_tag_count'] >= min_thresh]

    ko_only = ko_data[
        ko_data['naive_atac_tag_count'] < min_thresh]

    save_path = yzer.get_and_create_path(
        dirpath, 'Figures', 'Foxo1_group_overlaps')

    groups = [wt_only, both, ko_only]
    labels = ['WT only', 'WT and KO', 'Foxo1 KO only']

    if False:
        yzer.boxplot([gp['naive_foxo1_tag_count'] for gp in groups],
                     labels,
                     title='Foxo1 tags in ATAC-seq regions by group',
                     ylabel='Foxo1 peak tag count', save_dir=save_path,
                     show_plot=False)
        yzer.boxplot([gp['lcmv_d12_foxo1_tag_count'] for gp in groups],
                     labels,
                     title='LCMV d12 Foxo1 tags in ATAC-seq regions by group',
                     ylabel='Foxo1 peak tag count', save_dir=save_path,
                     show_plot=False)
        yzer.boxplot([gp['naive_h3k4me2_tag_count'] for gp in groups],
                     labels,
                     title='H3K4me2 tags in ATAC-seq regions by group',
                     ylabel='H3K4me2 region count', save_dir=save_path,
                     show_plot=False)
        yzer.boxplot([gp['lcmv_d12_h3k4me2_tag_count'] for gp in groups],
                     labels,
                     title='LCMV d12 H3K4me2 tags in ATAC-seq regions by group',
                     ylabel='H3K4me2 region count', save_dir=save_path,
                     show_plot=False)

    if False:
        for i, gp in enumerate(groups):
            yzer.piechart([sum(gp['naive_foxo1_tag_count'] >= min_thresh),
                           sum(gp['naive_foxo1_tag_count'] < min_thresh)],
                          ['With Foxo1', 'Without Foxo1'],
                          title='Co-occurrence with Foxo1- ' + labels[i],
                          save_dir=save_path, show_plot=False)
            yzer.piechart([sum(gp['lcmv_d12_foxo1_tag_count'] >= min_thresh),
                           sum(gp['lcmv_d12_foxo1_tag_count'] < min_thresh)],
                          ['With Foxo1', 'Without Foxo1'],
                          title='Co-occurrence with LCMV d12 Foxo1- ' +
                          labels[i],
                          save_dir=save_path, show_plot=False)
            yzer.histogram(gp['naive_foxo1_tag_count'].tolist(), bins=20,
                           title='Foxo1 peak tag count distribution- ' +
                           labels[i],
                           xlabel='Tag count in Foxo1 peak',
                           ylabel='Number of peaks',
                           save_dir=save_path, show_plot=False)

    if False:
        yzer.boxplot([gp[gp['naive_foxo1_tag_count'] > 0]['naive_foxo1_tag_count']
                      for gp in groups],
                     labels,
                     title='ATAC-seq regions with Foxo1 peaks by group',
                     ylabel='Foxo1 peak tag count', save_dir=save_path,
                     show_plot=False)

    # TCF1
    if True:
        for i, gp in enumerate(groups):
            yzer.piechart([sum(gp['tcf1_tag_count'] >= min_thresh),
                           sum(gp['tcf1_tag_count'] < min_thresh)],
                          ['With Tcf1', 'Without Tcf1'],
                          title='Co-occurrence with Tcf1- ' + labels[i],
                          save_dir=save_path, show_plot=False)
            yzer.piechart([sum((gp['tcf1_tag_count'] >= min_thresh) &
                               (gp['naive_foxo1_tag_count'] >= min_thresh)),
                           sum((gp['tcf1_tag_count'] < min_thresh) |
                               (gp['naive_foxo1_tag_count'] < min_thresh))],
                          ['With Tcf1 and Foxo1', 'Without Tcf1 or Foxo1'],
                          title='Co-occurrence with Tcf1 and Foxo1- ' +
                          labels[i],
                          save_dir=save_path, show_plot=True)
