'''
Created on Jan 14, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_sets')
    
    kla_col='kla_lfc'
   
    tss_only = False
    img_dirpath = yzer.get_and_create_path(dirpath, 'interaction_tag_count', tss_only and 'genic' or 'all_interactions')
    
    # File generated in novel_me2_sites
    enhancers = yzer.import_file(yzer.get_filename(data_dirpath,
                    'all_enhancers_with_me2_and_{0}interaction_stats.txt'.format(tss_only and 'tss_' or '')))
    transcripts = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_vectors.txt'))
    enhancers = enhancers.merge(transcripts, how='left', on='id')
    
    enhancers = enhancers[enhancers['total_interactions'] > 0]
    enhancers['dmso_tags_per_bp'] = enhancers['dmso_tag_count']/enhancers['length']
    enhancers['kla_tags_per_bp'] = enhancers['kla_tag_count']/enhancers['length']
    
    # Plot tag counts versus interactions.
    ax = yzer.scatterplot(enhancers, 
                     xcolname='dmso_tags_per_bp', ycolname='notx_interactions', log=True,  
                     title='Interactions in Notx as a function of GRO-seq tag counts in DMSO', 
                     xlabel='GRO-seq tags per bp in DMSO', 
                     ylabel='Number of interactions {0}in Notx'.format(tss_only and 'with gene TSSs ' or ''), 
                     show_2x_range=True, plot_regression=False, 
                     show_count=True, show_correlation=True, show_legend=False, 
                     save_dir=img_dirpath, show_plot=True)
    ax = yzer.scatterplot(enhancers, 
                     xcolname='kla_tags_per_bp', ycolname='kla_4h_interactions', log=True,  
                     title='Interactions in KLA 4h as a function of GRO-seq tag counts in KLA 1h', 
                     xlabel='GRO-seq tags per bp in KLA 1h', 
                     ylabel='Number of interactions {0}in KLA 4h'.format(tss_only and 'with gene TSSs ' or ''), 
                     show_2x_range=True, plot_regression=False, 
                     show_count=True, show_correlation=True, show_legend=False, 
                     save_dir=img_dirpath, show_plot=True)