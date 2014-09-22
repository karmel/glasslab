'''
Created on Sep 21, 2014

@author: karmel

The LCMV d12 H3K4me2 data has more background than the naive H3K4me2 data.
We want to try to normalize this, so we pretend that the total in-peak
tags for LCMV are the same fraction of the total data as the tags from 
the overlapping naive H3K4me2 peaks. First we figure out the fraction of the
tags-in-peaks that the overlapping naive peaks represent, 
and then scale the LCMV peaks such that both samples have comparable numbers of
tags-in-peaks relative to how many peaks there are.
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    yzer = SeqGrapher()

    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/' +\
        'Miscellaneous_Collaborations/Rodrigo_CD8s_2014_09/Normalization'
    dirpath = yzer.get_path(dirpath)

    lcmv_me2 = yzer.get_filename(dirpath, 'lcmv_h3k4me2_peaks.txt')
    lcmv_me2 = yzer.import_file(lcmv_me2)
    naive_me2 = yzer.get_filename(dirpath, 'naive_h3k4me2_peaks.txt')
    naive_me2 = yzer.import_file(naive_me2)

    total_lcmv_tags = sum(lcmv_me2['tag_count'].fillna(0))
    total_naive_tags = sum(naive_me2['tag_count'].fillna(0))
    overlapping_naive_tags = sum(lcmv_me2['naive_h3k4me2_tag_count'].fillna(0))
    naive_fraction = overlapping_naive_tags / total_naive_tags
    print(naive_fraction)
    lcmv_me2['tag_count'] = lcmv_me2['tag_count'] * naive_fraction * \
        total_naive_tags / total_lcmv_tags

    output_file = yzer.get_filename(dirpath, 'lcmv_h3k4me2_peaks_norm.txt')
    lcmv_me2.to_csv(output_file, sep='\t', header=True, index=False)
