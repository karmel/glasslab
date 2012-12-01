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
        img_dirpath = yzer.get_and_create_path(dirpath, 'fold_change_per_domain', 'up_in_kla',
                                               'rep{0}'.format(rep))
        
        kla_key = 'kla_{0}_lfc'.format(rep)
        dex_kla_key = 'dex_over_kla_{0}_lfc'.format(rep)
        
        data = data[data[kla_key] > 1]
        
        shuffled = data['domain_id'].values.copy()
        shuffle(shuffled)
        data['shuffled_domain_id'] = shuffled 
        
        data['up_in_kla'] = data[kla_key] > 1
        data['repressed'] = data[dex_kla_key] <= -.58
        data['transrepressed'] = (data[kla_key] > 1) & (data[dex_kla_key] <= -.58)
        data['count'] = ~data[kla_key].isnull()
        
        grouped = data.groupby(by='domain_id', as_index=False).sum()
        shuffled_grouped = data.groupby(by='shuffled_domain_id', as_index=False).sum()
        
        grouped = grouped[grouped['domain_id'] != 0]
        shuffled_grouped = shuffled_grouped[shuffled_grouped['shuffled_domain_id'] != 0]
        '''
        grouped['kla_ratio'] = grouped['up_in_kla']/grouped['count']
        grouped = grouped.sort(['kla_ratio']).reset_index(drop=True)
        grouped['idx'] = grouped.index
        
        shuffled_grouped['kla_ratio'] = shuffled_grouped['up_in_kla']/shuffled_grouped['count']
        shuffled_grouped = shuffled_grouped.sort(['kla_ratio']).reset_index(drop=True)
        shuffled_grouped['idx'] = shuffled_grouped.index
        
        ax = yzer.scatterplot(shuffled_grouped, 'idx', 'kla_ratio', 
                         color='green', label='Shuffled Data'.format(rep),
                         show_2x_range=False, plot_regression=False, 
                         show_count=False, show_correlation=False, 
                         show_legend=False, show_plot=False)
        
        ax = yzer.scatterplot(grouped, 'idx', 'kla_ratio', 
                         title='Up in KLA {0} Percentage by HiC Domain'.format(rep), 
                         xlabel='Ordered Index', ylabel='Percent of transcripts up in KLA', 
                         color='blue', label='Replicate {0} Data'.format(rep),
                         show_2x_range=False, plot_regression=False, 
                         show_count=False, show_correlation=False, 
                         show_legend=True, save_dir=img_dirpath, show_plot=True, ax=ax)
        
        '''
        grouped['repressed_ratio'] = grouped['repressed']/grouped['count']
        grouped = grouped.sort(['repressed_ratio']).reset_index(drop=True)
        grouped['idx'] = grouped.index
        
        shuffled_grouped['repressed_ratio'] = shuffled_grouped['repressed']/shuffled_grouped['count']
        shuffled_grouped = shuffled_grouped.sort(['repressed_ratio']).reset_index(drop=True)
        shuffled_grouped['idx'] = shuffled_grouped.index
        
        ax = yzer.scatterplot(shuffled_grouped, 'idx', 'repressed_ratio', 
                         color='green', label='Shuffled Data'.format(rep),
                         show_2x_range=False, plot_regression=False, 
                         show_count=False, show_correlation=False, 
                         show_legend=False, show_plot=False)
        
        ax = yzer.scatterplot(grouped, 'idx', 'repressed_ratio', 
                         title='Repressed by Dex {0} Percentage by HiC Domain'.format(rep), 
                         xlabel='Ordered Index', ylabel='Percent of transcripts down in Dex', 
                         color='blue', label='Replicate {0} Data'.format(rep),
                         show_2x_range=False, plot_regression=False, 
                         show_count=False, show_correlation=False, 
                         show_legend=True, save_dir=img_dirpath, show_plot=True, ax=ax)
       
       
        grouped['transrepressed_ratio'] = grouped['transrepressed']/grouped['count']
        grouped = grouped.sort(['transrepressed_ratio']).reset_index(drop=True)
        grouped['idx'] = grouped.index
        
        shuffled_grouped['transrepressed_ratio'] = shuffled_grouped['transrepressed']/shuffled_grouped['count']
        shuffled_grouped = shuffled_grouped.sort(['transrepressed_ratio']).reset_index(drop=True)
        shuffled_grouped['idx'] = shuffled_grouped.index
        
        ax = yzer.scatterplot(shuffled_grouped, 'idx', 'transrepressed_ratio', 
                         color='green', label='Shuffled Data'.format(rep),
                         show_2x_range=False, plot_regression=False, 
                         show_count=False, show_correlation=False, 
                         show_legend=False, show_plot=False)
        
        ax = yzer.scatterplot(grouped, 'idx', 'transrepressed_ratio', 
                         title='Transepressed by Dex {0} Percentage by HiC Domain'.format(rep), 
                         xlabel='Ordered Index', ylabel='Percent of transcripts transrepressed in Dex', 
                         color='blue', label='Replicate {0} Data'.format(rep),
                         show_2x_range=False, plot_regression=False, 
                         show_count=False, show_correlation=False, 
                         show_legend=True, save_dir=img_dirpath, show_plot=True, ax=ax)
       
       