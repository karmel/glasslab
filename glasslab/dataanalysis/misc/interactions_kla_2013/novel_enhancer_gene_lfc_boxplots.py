'''
Created on Jan 9, 2013

@author: karmel

Do novel interactions in KLA correlate with increased KLA expression?
Are lost interactions responsible for the skew downwards?
'''

from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_sets')
    img_dirpath = yzer.get_and_create_path(dirpath, 'novel_enhancer_gene_lfc')
    
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
    
    interactions['hash'] = interactions.apply(lambda row: '{0}.{1}'.format(row['id'],row['id_2']), axis=1)
    
    pairs = {}
    pairs['notx'] = interactions[interactions['sequencing_run_id'] == 765]
    pairs['kla_30m'] = interactions[interactions['sequencing_run_id'] == 766]
    pairs['kla_4h'] = interactions[interactions['sequencing_run_id'] == 773]
    
    
    # Get exclusive sets of notx versus KLA-- note that it is okay for pairs to be shared in diff KLA
    pairs['only_notx'] = pairs['notx'][(~pairs['notx']['hash'].isin(pairs['kla_30m']['hash']))
                                        & (~pairs['notx']['hash'].isin(pairs['kla_4h']['hash']))]
    pairs['only_kla_30m'] = pairs['kla_30m'][(~pairs['kla_30m']['hash'].isin(pairs['notx']['hash']))]
    pairs['only_kla_4h'] = pairs['kla_4h'][(~pairs['kla_4h']['hash'].isin(pairs['notx']['hash']))]
    
    # Draw boxplots for KLA LFC for enhancers and genes in each set of pairs.
    labels = ['Enhancers, notx', 'Genes, notx',
              'Enhancers, KLA 30m', 'Genes, KLA 30m',
              'Enhancers, KLA 4h', 'Genes, KLA 4h',]
    vals = [pairs['only_notx'][kla_col + '_2'], pairs['only_notx'][kla_col],
            pairs['only_kla_30m'][kla_col + '_2'], pairs['only_kla_30m'][kla_col],
            pairs['only_kla_4h'][kla_col + '_2'], pairs['only_kla_4h'][kla_col],
            ]
    
    labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(labels, vals)]
    title = 'KLA 1h Log Fold Change for Interactions Unique to notx or KLA (Score >= 10)'
    ax = yzer.boxplot(vals, labels, 
                     title=title, xlabel='Transcript subset, sample of interaction', 
                     ylabel='log2(KLA 1h GRO-seq/notx GRO-seq)', 
                     show_outliers=False, show_plot=True, wide=True,
                     save_dir=img_dirpath)