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
    img_dirpath = yzer.get_and_create_path(dirpath, 'novel_interactions_kla_lfc','all_interactions')
    
    interactions = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_enhancer_with_anything_with_me2_inc_me2_counts.txt'))
    interactions = interactions[interactions['count'] > 1]

    all_transcripts = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_vectors.txt'))
    
    kla_col='kla_lfc'
    
    transcripts = all_transcripts[['id', kla_col]]
        
    # Associate gene id
    interactions = interactions.merge(transcripts, how='left', on='id')
    
    transcripts['id_2'] = transcripts['id']
    transcripts = transcripts.drop(['id'], axis=1)
    interactions = interactions.merge(transcripts, how='left', on='id_2', suffixes=['','_2'])
    
    interactions = interactions.fillna(0)
    
    # Key on peak id, not enhancer id, which could be bidirectional
    #interactions['id_2'] = interactions['h3k4me2_id']
    interactions['hash'] = interactions.apply(lambda row: '{0}.{1}'.format(row['id'],row['id_2']), axis=1)
    
    
    
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
              'Partner; newly interacting', 'Partner; switched interaction', 'Partner; existing',]
    vals = [pairs['only_kla_30m'][kla_col + '_2'], pairs['only_kla_30m_ixn'][kla_col + '_2'], pairs['shared_kla_30m'][kla_col + '_2'],  
             pairs['only_kla_30m'][kla_col], pairs['only_kla_30m_ixn'][kla_col], pairs['shared_kla_30m'][kla_col],
            ]
    
    labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(labels, vals)]
    title = 'KLA 1h LFC for Novel Interactions in KLA 30m'
    ax = yzer.boxplot(vals, labels, 
                     title=title, xlabel='Transcript subset; interaction state in KLA 30m versus notx', 
                     ylabel='log2(KLA 1h GRO-seq/notx GRO-seq)', 
                     show_outliers=False, show_plot=True, wide=True,
                     save_dir=img_dirpath)