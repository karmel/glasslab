'''
Created on Mar 11, 2013

@author: karmel

Note: Made font.weight = bold and axes.titlesize = 24, font.size = 16 in matplotlibrc
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/NAR_review_data/vs_homer'
    dirpath = yzer.get_path(dirpath)
    
    img_dirpath = yzer.get_and_create_path(dirpath, 'scatterplots')
    
    data = yzer.import_file(yzer.get_filename(dirpath,'tag_count_by_refseq.txt'))
    data['sum'] = nonzero(data['sum'].fillna(0))
    
    homer_data = yzer.import_file(yzer.get_filename(dirpath,'RNA_GroSeq_CountsGenes.txt'))
    homer_data['sequence_identifier'] = homer_data['Gene ID']
    homer_data['homer_tag_count'] = nonzero(homer_data['ThioMac-GroSeq-notx-110513/ genes (Total: 12166480.0) normFactor 0.82'].fillna(0))
    homer_data = homer_data[['sequence_identifier','homer_tag_count']]
    
    merged = data.merge(homer_data, how='inner',on='sequence_identifier')
    merged = merged.fillna(1)
    
    if True:
        ax = yzer.scatterplot(merged, 
                     xcolname='homer_tag_count', ycolname='sum', log=True,  
                     title='RefSeq Tag Count via Homer and Vespucci', 
                     xlabel='Tag Count in Homer',
                     ylabel='Tag Count in Vespucci',  
                     show_2x_range=True, plot_regression=False, set_limits=True,
                     show_count=True, show_correlation=True, show_legend=False, 
                     save_dir=img_dirpath, show_plot=True)

    merged['ratio'] = merged['sum']/merged['homer_tag_count']
    merged = merged.sort('ratio')
    print merged.head(10)
    print merged.tail(20)