'''
Created on Jan 3, 2013

@author: karmel

Plot gen-enhancer me2 LFC; do we see correlation?
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero
import numpy

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_sets')
    img_dirpath = yzer.get_and_create_path(dirpath, 'gene_enhancer_me2_lfc','scatterplots')
    
    interactions = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_enhancer_with_anything_with_me2_inc_me2_counts.txt'))
    interactions = interactions[interactions['count'] > 1]
    all_transcripts = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_vectors.txt'))
    
    for me2_timepoint in ('6h','24h'):
        me2_col = 'me2_{0}_ratio'.format(me2_timepoint)
        kla_col = 'kla_lfc'
        col_set = [me2_col + '_2', kla_col + '_2', kla_col, me2_col]
        
        interactions[me2_col] = numpy.log2(nonzero(interactions['me2_kla_{0}_tag_count'.format(me2_timepoint)])/\
                                            nonzero(interactions['me2_notx_tag_count']))
        interactions[me2_col + '_2'] = numpy.log2(nonzero(interactions['me2_kla_{0}_tag_count_2'.format(me2_timepoint)])/\
                                            nonzero(interactions['me2_notx_tag_count_2']))
        
        
        transcripts = all_transcripts[['id', kla_col]]
        
        # Associate gene id
        interactions = interactions.merge(transcripts, how='left', on='id')
        
        transcripts['id_2'] = transcripts['id']
        transcripts = transcripts.drop(['id'], axis=1)
        interactions = interactions.merge(transcripts, how='left', on='id_2', suffixes=['','_2'])
        
        interactions = interactions.fillna(0)
        
        pairs = {}
        pairs['notx'] = interactions[interactions['sequencing_run_id'] == 765][col_set]
        pairs['kla_30m'] = interactions[interactions['sequencing_run_id'] == 766][col_set]
        pairs['kla_4h'] = interactions[interactions['sequencing_run_id'] == 773][col_set]
    
        for key, val_set  in pairs.iteritems():
            
            val_set = val_set[(val_set[me2_col + '_2'].abs() >= .58)]
            print key, len(val_set)
            print val_set.corr()
            
            if True:
                val_set = val_set.sort(columns=me2_col + '_2', ascending=True)
                change_subdir = yzer.get_and_create_path(dirpath, 'gene_enhancer_me2_lfc','changing_me2')
                
                # Fake col for .cdt file
                # We want 'id    weight    enhancer_lfc    gene_lfc'
                val_set['weight'] = 1.0
                f = open(yzer.get_filename(change_subdir, '{0}_for_{1}_pairs.cdt'.format(me2_col, key)),'w')
                val_set.to_csv(f, sep='\t', header=False, index=True, cols=(['weight'] + col_set))
                
                if False:
                    ax = yzer.scatterplot(val_set, 
                         xcolname=me2_col + '_2', ycolname=me2_col, log=False, color='blue', 
                         title='Log fold change of genes and interacting enhancers in {0}: {1}, enhancer 1.5x fold changed'.format(
                                                                        key.replace('_',' '), kla_col), 
                         xlabel='Enhancer LFC', ylabel='Gene LFC', 
                         show_2x_range=False, plot_regression=True, 
                         show_count=True, show_correlation=True, show_legend=False, 
                         save_dir=img_dirpath, show_plot=False)