'''
Created on Jan 9, 2013

@author: karmel

Do novel interactions gain or lose me2? 
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_sets')
    img_dirpath = yzer.get_and_create_path(dirpath, 'novel_enhancer_me2_change')
    
    interactions = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_refseq_with_me2_inc_me2_counts.txt'))
    interactions = interactions[interactions['count'] > 1]

    
    interactions = interactions.fillna(0)
    
    interactions['hash'] = interactions.apply(lambda row: '{0}.{1}'.format(row['id'],row['id_2']), axis=1)
    
    for me2_timepoint in ('6h', '24h'):
        interactions['me2_ratio'] = nonzero(interactions['me2_kla_{0}_tag_count'.format(me2_timepoint)])/\
                                            nonzero(interactions['me2_notx_tag_count'])
        interactions['me2_ratio_2'] = nonzero(interactions['me2_kla_{0}_tag_count_2'.format(me2_timepoint)])/\
                                            nonzero(interactions['me2_notx_tag_count_2'])
        
        
        col = 'me2_ratio'
        pairs = {}
        pairs['notx'] = interactions[interactions['sequencing_run_id'] == 765]
        pairs['kla_30m'] = interactions[interactions['sequencing_run_id'] == 766]
        pairs['kla_4h'] = interactions[interactions['sequencing_run_id'] == 773]
        
        
        # Enhancer is totally new interactor in 4h KLA
        pairs['only_notx'] = pairs['notx'][(~pairs['notx']['id_2'].isin(pairs['kla_30m']['id_2']))
                                            & (~pairs['notx']['id_2'].isin(pairs['kla_4h']['id_2']))]
        pairs['only_kla_30m'] = pairs['kla_30m'][(~pairs['kla_30m']['id_2'].isin(pairs['notx']['id_2']))]
        pairs['only_kla_4h'] = pairs['kla_4h'][(~pairs['kla_4h']['id_2'].isin(pairs['notx']['id_2']))]
        # Interaction is switched in 4h KLA
        pairs['only_notx_ixn'] = pairs['notx'][(~pairs['notx']['hash'].isin(pairs['kla_30m']['hash']))
                                            & (~pairs['notx']['hash'].isin(pairs['kla_4h']['hash']))
                                            & (~pairs['notx']['hash'].isin(pairs['only_notx']['hash']))]
        pairs['only_kla_30m_ixn'] = pairs['kla_30m'][(~pairs['kla_30m']['hash'].isin(pairs['notx']['hash']))
                                             & (~pairs['kla_30m']['hash'].isin(pairs['only_kla_30m']['hash']))]
        pairs['only_kla_4h_ixn'] = pairs['kla_4h'][(~pairs['kla_4h']['hash'].isin(pairs['notx']['hash']))
                                                & (~pairs['kla_4h']['hash'].isin(pairs['only_kla_4h']['hash']))]
        
        # Else, existing interaction
        pairs['shared_notx'] = pairs['notx'][(pairs['notx']['hash'].isin(pairs['kla_30m']['hash']))
                                            | (pairs['notx']['hash'].isin(pairs['kla_4h']['hash']))]
        pairs['shared_kla_30m'] = pairs['kla_30m'][(pairs['kla_30m']['hash'].isin(pairs['notx']['hash']))]
        pairs['shared_kla_4h'] = pairs['kla_4h'][(pairs['kla_4h']['hash'].isin(pairs['notx']['hash']))]
        
        # Draw boxplots for KLA LFC for enhancers and genes in each set of pairs.
        labels = ['Enhancers; newly interacting', 'Enhancers; switched interaction', 'Enhancers; existing', 
                  'Genes; newly interacting', 'Genes; switched interaction', 'Genes; existing',]
        vals = [pairs['only_kla_4h'][col + '_2'], pairs['only_kla_4h_ixn'][col + '_2'], pairs['shared_kla_4h'][col + '_2'],  
                 pairs['only_kla_4h'][col], pairs['only_kla_4h_ixn'][col], pairs['shared_kla_4h'][col],
                ]
        
        labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(labels, vals)]
        title = 'KLA {0} H3K4me2 Tag Count Ratio for Novel Interactions in KLA 4h'.format(me2_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Transcript subset; interaction state in KLA 4h versus notx', 
                         ylabel='H3K4me2 Ratio: me2 in KLA {0}/me2 in notx'.format(me2_timepoint), 
                         show_outliers=False, show_plot=True, wide=True,
                         save_dir=img_dirpath)