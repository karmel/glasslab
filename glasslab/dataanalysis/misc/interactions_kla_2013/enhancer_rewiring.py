'''
Created on Jan 30, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from collections import OrderedDict

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_sets')
    img_dirpath = yzer.get_and_create_path(dirpath, 'enhancer_rewiring')
    
    interactions = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_refseq_with_me2.txt'))
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
    
    
    subsets = OrderedDict((('All', interactions.copy()),
                          ('Up-reg genes', interactions[interactions['kla_lfc'] > 1]),
                          ('Down-reg genes', interactions[interactions['kla_lfc'] < -1]),
                          ('Up-reg enh', interactions[interactions['kla_lfc_2'] > 1]),
                          ('Down-reg enh', interactions[interactions['kla_lfc_2'] < -1]),
                          )) 
    
    for key, sub in subsets.iteritems():
        print '\nCounts for ' + key
        print '============================='
        pairs = {}
        pairs['notx'] = sub[sub['sequencing_run_id'] == 765]
        pairs['kla_30m'] = sub[sub['sequencing_run_id'] == 766]
        pairs['kla_4h'] = sub[sub['sequencing_run_id'] == 773]
        timepoints = pairs.keys()
        
        # Interaction is switched in 4h KLA
        pairs['only_notx'] = pairs['notx'][(~pairs['notx']['hash'].isin(pairs['kla_4h']['hash']))
                                           & (~pairs['notx']['hash'].isin(pairs['kla_30m']['hash']))]
        pairs['only_kla_30m'] = pairs['kla_30m'][(~pairs['kla_30m']['hash'].isin(pairs['notx']['hash']))]
        pairs['only_kla_4h'] = pairs['kla_4h'][(~pairs['kla_4h']['hash'].isin(pairs['notx']['hash']))]
        
        # Else, existing interaction
        pairs['shared_notx'] = pairs['notx'][(pairs['notx']['hash'].isin(pairs['kla_30m']['hash']))
                                            | (pairs['notx']['hash'].isin(pairs['kla_4h']['hash']))]
        pairs['shared_kla_30m'] = pairs['kla_30m'][(pairs['kla_30m']['hash'].isin(pairs['notx']['hash']))]
        pairs['shared_kla_4h'] = pairs['kla_4h'][(pairs['kla_4h']['hash'].isin(pairs['notx']['hash']))]
    
        for t in timepoints:
            print '\n........................\nCounts for ' + t
            print 'Number of interactions (total, unique, shared):'
            print len(pairs[t]), len(pairs['only_' + t]), len(pairs['shared_' + t])
            
            print 'Number of implicated genes (total, unique, shared):'
            print len(pairs[t]['id'].unique()), len(pairs['only_' + t]['id'].unique()), len(pairs['shared_' + t]['id'].unique())
            
            
            
        