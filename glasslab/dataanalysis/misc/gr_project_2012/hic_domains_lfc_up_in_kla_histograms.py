'''
Created on Nov 26, 2012

@author: karmel


'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from random import shuffle


if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/hic_domains'
    dirpath = yzer.get_path(dirpath)
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
    data = data.fillna(0)
    
    for rep in (4,3,1):
        img_dirpath = yzer.get_and_create_path(dirpath, 'lfc_histograms',
                                               'rep{0}'.format(rep))
        
        kla_key = 'kla_{0}_lfc'.format(rep)
        dex_kla_key = 'dex_over_kla_{0}_lfc'.format(rep)
        
        shuffled = data['domain_id'].values.copy()
        shuffle(shuffled)
        data['shuffled_domain_id'] = shuffled 
        
        data = data[data[kla_key] >= 1]
        
        grouped = data.groupby(by='domain_id', as_index=False).mean()
        shuffled_grouped = data.groupby(by='shuffled_domain_id', as_index=False).mean()
        
        grouped = grouped[grouped['domain_id'] != 0]
        shuffled_grouped = shuffled_grouped[shuffled_grouped['shuffled_domain_id'] != 0]
        
        
        
        ax = yzer.histogram(grouped[dex_kla_key], bins=50,
                         label='Replicate {0} Data'.format(rep),
                         show_legend=False, show_plot=False)
        
        ax = yzer.histogram(shuffled_grouped[dex_kla_key], bins=50,
                         title='Dex+KLA over KLA LFC {0} by HiC Domain'.format(rep), 
                         xlabel='Mean Dex+KLA log fold change for transcripts up in KLA', ylabel='Count of Domains', 
                         color='black', fill=False,label='Shuffled Data'.format(rep),
                         show_legend=True, save_dir=img_dirpath, show_plot=False, ax=ax)
        
