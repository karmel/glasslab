'''
Created on Jan 3, 2013

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/gene_enhancer_lfc'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'scatterplots')
    
    interactions = yzer.import_file(yzer.get_filename(dirpath,'refseq_with_distal_norm_count_gte_2.txt'))
    all_transcripts = yzer.import_file(yzer.get_filename(dirpath,'transcript_vectors.txt'))
    
    for x in xrange(0,5):
        if not x: kla_col = 'kla_lfc'
        else: kla_col = 'kla_{0}_lfc'.format(x)
        
        transcripts = all_transcripts[['id', kla_col]]
        
        # Associate gene id
        interactions = interactions.merge(transcripts, how='left', on='id')
        
        transcripts['id_2'] = transcripts['id']
        transcripts = transcripts.drop(['id'], axis=1)
        interactions = interactions.merge(transcripts, how='left', on='id_2', suffixes=['','_2'])
        
        interactions = interactions.fillna(0)
        
        pairs = {}
        pairs['notx'] = interactions[interactions['sequencing_run_id'] == 765][[kla_col + '_2', kla_col]]
        pairs['kla_30m'] = interactions[interactions['sequencing_run_id'] == 766][[kla_col + '_2', kla_col]]
        pairs['kla_4h'] = interactions[interactions['sequencing_run_id'] == 773][[kla_col + '_2', kla_col]]
    
        for key, val_set  in pairs.iteritems():
            
            val_set = val_set[(val_set[kla_col + '_2'].abs() >= .58)]
            print key, len(val_set)
            print val_set.corr()
            
            val_set = val_set.sort(columns=kla_col + '_2', ascending=True)
            f = open(yzer.get_filename(dirpath, 'changing_only_enh', '{0}_for_{1}_pairs.txt'.format(kla_col, key)),'w')
            val_set.to_csv(f, sep='\t', header=False, index=True)
            
            ax = yzer.scatterplot(val_set, 
                 xcolname=kla_col + '_2', ycolname=kla_col, log=False, color='blue', 
                 title='Log fold change of genes and interacting enhancers in {0}: {1}, enhancer 1.5x fold changed'.format(
                                                                key.replace('_',' '), kla_col), 
                 xlabel='Enhancer LFC', ylabel='Gene LFC', 
                 show_2x_range=False, plot_regression=True, 
                 show_count=True, show_correlation=True, show_legend=False, 
                 save_dir=img_dirpath, show_plot=False)