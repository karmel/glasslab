'''
Created on Sep 20, 2014

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    yzer = SeqGrapher()

    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/' +\
        'Miscellaneous_Collaborations/Rodrigo_CD8s_2014_09/Enhancers'
    dirpath = yzer.get_path(dirpath)

    if True:
        atac = yzer.get_filename(
            dirpath, 'naive_atac', 'naive_atac_enhancers.txt')
        me2 = yzer.get_filename(
            dirpath, 'naive_h3k4me2', 'naive_h3k4me2_enhancers.txt')
        ac = yzer.get_filename(
            dirpath, 'naive_h3k27ac', 'naive_h3k27ac_enhancers.txt')

        atac = yzer.import_file(atac).fillna(0)
        me2 = yzer.import_file(me2).fillna(0)

        peak_thresh, region_thresh = 10, 20
        atac = atac[atac['tag_count'] >= peak_thresh]
        me2 = me2[me2['tag_count'] >= region_thresh]

        # Make sure to count each separately, as many ATAC peaks
        # can be subsumed by a single H3K4me2 peak
        atac_only = atac[(atac['naive_h3k4me2_tag_count'] < region_thresh)]
        atac_me2 = atac[(atac['naive_h3k4me2_tag_count'] >= region_thresh)]
        me2_only = me2[(me2['naive_atac_tag_count'] < peak_thresh)]
        me2_atac = me2[(me2['naive_atac_tag_count'] >= peak_thresh)]

        print('ATAC only: ', len(atac_only))
        print('ATAC with H3K4me2: ', len(atac_me2))
        print('H3K4me2 only: ', len(me2_only))
        print('H3K4me2 with ATAC: ', len(me2_atac))

        yzer.piechart([len(atac_only), len(atac_me2)],
                      ['ATAC only', 'ATAC with H3K4me2'],
                      title='ATAC-seq region overlaps',
                      save_dir=yzer.get_and_create_path(dirpath, 'Figures'))

        yzer.piechart([len(me2_only), len(me2_atac)],
                      ['H3K4me2 only', 'H3K4me2 with ATAC'],
                      title='H3K4me2 overlaps',
                      save_dir=yzer.get_and_create_path(dirpath, 'Figures'))

        yzer.boxplot([atac_only['tag_count'], atac_me2['tag_count']],
                     ['ATAC only', 'ATAC with H3K4me2'],
                     title='ATAC-seq tag counts by H3K4me2 overlap',
                     xlabel='Group', ylabel='Peak tag count',
                     save_dir=yzer.get_and_create_path(dirpath, 'Figures'))
        yzer.boxplot([me2_only['tag_count'], me2_atac['tag_count']],
                     ['H3K4me2 only', 'H3K4me2 with ATAC'],
                     title='H3K4me2 tag counts by ATAC-seq overlap',
                     xlabel='Group', ylabel='Peak tag count',
                     save_dir=yzer.get_and_create_path(dirpath, 'Figures'))
        yzer.histogram(atac_only['tag_count'].tolist(), bins=20,
                       title='ATAC-seq-only peak tag count distribution',
                       xlabel='Tag count in peak', ylabel='Number of peaks',
                       save_dir=yzer.get_and_create_path(dirpath, 'Figures'))
        yzer.histogram(me2_only['tag_count'].tolist(), bins=20,
                       title='H3K4me2-only peak tag count distribution',
                       xlabel='Tag count in peak', ylabel='Number of peaks',
                       save_dir=yzer.get_and_create_path(dirpath, 'Figures'))
