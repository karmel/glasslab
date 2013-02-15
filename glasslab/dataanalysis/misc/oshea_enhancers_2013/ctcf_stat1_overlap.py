'''
Created on Feb 14, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/Oshea_enhancers/ctcf_stat1_overlap'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'figures')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'ctcf_with_stat1_binding.txt')).fillna(0)
    with_stat1 = data[data['p2_tag_count'] > 0]
    without_stat1 = data[data['p2_tag_count'] == 0]

    if False:
        ax = yzer.piechart([len(with_stat1), len(without_stat1)], 
                      ['CTCF sites with STAT1', 'CTCF sites without STAT1'], 
                      title='DP Thymocyte CTCF Sites with STAT1 in Th1 Cells', 
                      save_dir=img_dirpath, 
                      show_plot=True)
    data['tag_count_nonzero'] = nonzero(data['tag_count'])
    data['p2_tag_count_nonzero'] = nonzero(data['p2_tag_count'])
    ax = yzer.scatterplot(data, 'tag_count_nonzero', 'p2_tag_count_nonzero',
                            xlabel='CTCF Tag Count', ylabel='Stat1 Tag Count',
                            log=True, color='blue', 
                            title='Tags in CTCF Peaks versus Overlapping Stat1 Peaks',
                            show_2x_range=False, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)