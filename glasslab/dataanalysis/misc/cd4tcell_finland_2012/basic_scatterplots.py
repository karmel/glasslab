'''
Created on Nov 7, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.dataanalysis.misc.cd4tcell_finland_2012.resources import comparison_sets,\
    pretty_names

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells_Finland_2012/Analysis_2013_02'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'with_me3', 'basic_scatterplots')

    data = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    data = data.fillna(0)
    data = data[data['naive_me3_tag_count'] + data['act_me3_tag_count'] > 0]
                            
    for key1, key2, norm_factor in comparison_sets:
        name1 = pretty_names[key1[:-1]] + key1[-1:]
        name2 = pretty_names[key2[:-1]] + key2[-1:]
        
        data_normed = yzer.normalize(data, key2 + '_tag_count', norm_factor)
        ax = yzer.scatterplot(data_normed, key1 + '_tag_count', key2 + '_tag_count_norm',
                    log=True, color='blue', 
                    title='{0} versus {1} Normalized Tag Counts'.format(name1, name2),
                    xlabel='{0} tags in RefSeq transcripts'.format(name1), 
                    ylabel='{0} tags in RefSeq transcripts, normalized'.format(name2),
                    add_noise=False,
                    show_2x_range=True, show_legend=False, 
                    plot_regression=False,
                    show_count=True, show_correlation=True, 
                    save_dir=img_dirpath, show_plot=False)