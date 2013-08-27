'''
Created on Aug 25, 2013

@author: karmel

Plot supplementary figure showing Hah et al error rates 
against MAX EDGE values when Vespucci is built without knowledge
of RefSeq boundaries.
 
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/NAR_review_data/no_refseq'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'plots')

    ax = yzer.set_up_plot()
    title = 'Benchmarking without Foreknowledge of RefSeq'
    yzer.add_title(title, ax)
    yzer.add_axis_labels('MAX_EDGE value', 
                         'Error rate defined by Hah et al. (%)')
    
    max_edges = [100, 500, 1000, 4000, 5000, 5500, 6000, 10000]
    error_rates = [0.388551822833, 0.372390444765, 
                   0.263807982126, 0.124663089396, 
                   0.121784970634, 0.121807917409, 
                   0.123263849815, 0.142530838464]
    error_pcts = [e*100 for e in error_rates]
    yzer.plot(max_edges, error_pcts, '-o')
    yzer.save_plot_with_dir(save_dir=img_dirpath, title=title)
    yzer.show_plot()
    