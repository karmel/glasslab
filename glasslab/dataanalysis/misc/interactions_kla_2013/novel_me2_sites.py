'''
Created on Jan 11, 2013

@author: karmel

What do enhancers that are gaining methyl with KLA look like?

'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero
from collections import OrderedDict

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_sets')
    
    kla_col='kla_lfc'
   
    tss_only = False
    img_dirpath = yzer.get_and_create_path(dirpath, 'novel_me2_sites', tss_only and 'genic' or 'all_interactions','ratio_4')
 
    if False:
        enhancers = yzer.import_file(yzer.get_filename(data_dirpath,'all_distal_enhancers_inc_me2.txt'))
        
        all_transcripts = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_vectors.txt'))
        transcripts = all_transcripts[['id', kla_col]]
        enhancers = enhancers.merge(transcripts, how='left', on='id')
        
        if tss_only: interactions = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_refseq_with_me2.txt'))
        else: interactions = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_enhancer_with_anything.txt'))
        interactions = interactions[interactions['count'] > 1]
        interactions = interactions.fillna(0)
        interactions['hash'] = interactions.apply(lambda row: '{0}.{1}'.format(row['id'],row['id_2']), axis=1)
        
        # Initialize stats columns
        for key in ('total_interactions','notx_interactions',
                    'kla_30m_interactions','kla_4h_interactions',
                    'kla_30m_over_notx','kla_4h_over_notx'): 
            enhancers[key] = 0.0 # Don't use int here or else pandas forces the int type.
        
        # Add some stats
        for index, row in enhancers.iterrows():
            relevant_intxns = interactions[interactions['id_2'] == row['id']]
            # Straight counts
            enhancers.ix[index,'total_interactions'] = sum(relevant_intxns.get('norm_count',[]))
            enhancers.ix[index,'notx_interactions'] = sum(relevant_intxns[relevant_intxns['sequencing_run_id'] == 765].get('norm_count',[]))
            enhancers.ix[index,'kla_30m_interactions'] = sum(relevant_intxns[relevant_intxns['sequencing_run_id'] == 766].get('norm_count',[]))
            enhancers.ix[index,'kla_4h_interactions'] = sum(relevant_intxns[relevant_intxns['sequencing_run_id'] == 773].get('norm_count',[]))
            enhancers.ix[index,'kla_30m_over_notx'] =  (enhancers.ix[index,'kla_30m_interactions'] or .1)/\
                                                    (enhancers.ix[index,'notx_interactions'] or .1)
            enhancers.ix[index,'kla_4h_over_notx'] =  (enhancers.ix[index,'kla_4h_interactions'] or .1)/\
                                                    (enhancers.ix[index,'notx_interactions'] or .1)
            
        
        # Output so that we don't have to recompute that every time.
        enhancers.to_csv(yzer.get_filename(data_dirpath, 'all_enhancers_with_me2_and_interaction_stats.txt'), 
                         sep='\t', header=True, index=False)
    
    enhancers = yzer.import_file(yzer.get_filename(data_dirpath,
                    'all_enhancers_with_me2_and_{0}interaction_stats.txt'.format(tss_only and 'tss_' or '')))
    
    col = 'me2_ratio'
    for me2_timepoint in ('6h', '24h'):
        enhancers[col] = nonzero(enhancers['me2_kla_{0}_tag_count_2'.format(me2_timepoint)])/\
                                    nonzero(enhancers['me2_notx_tag_count_2'])
        
            
        sets = OrderedDict()
        sets['2x me2 in KLA {0}'.format(me2_timepoint)] = enhancers[enhancers[col] > 4]
        sets['No change me2 in KLA {0}'.format(me2_timepoint)] = enhancers[(enhancers[col] >= .5) & (enhancers[col] <= 2)]
        sets['1/2 me2 in KLA {0}'.format(me2_timepoint)] = enhancers[enhancers[col] < .25]
        
        labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(sets.keys(), sets.values())]
        
        vals = [v[kla_col] for v in sets.values()]
        title = 'KLA 1h LFC by Ratio of KLA {0} to Notx H3K4me2 for All Enhancers'.format(me2_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='All distal enhancers by relative amount of H3K4me2', 
                         ylabel='log2(KLA 1h GRO-seq/notx GRO-seq)', 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        
        # Now limit to interacting enhancers
        for key, val in sets.iteritems(): sets[key] = val[val['total_interactions'] > 0]
        
        labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(sets.keys(), sets.values())]
        
        vals = [v[kla_col] for v in sets.values()]
        title = 'KLA 1h Log Fold Change by Ratio of KLA {0} to Notx H3K4me2'.format(me2_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by relative amount of H3K4me2', 
                         ylabel='log2(KLA 1h GRO-seq/notx GRO-seq)', 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        
        vals = [v['total_interactions'] for v in sets.values()]
        title = 'Total Interactions with Genes by Ratio of KLA {0} to Notx H3K4me2'.format(me2_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by relative amount of H3K4me2', 
                         ylabel='Count of interactions {0}in all conditions'.format(tss_only and 'with gene TSSs ' or ''), 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        
        
        vals = [v['notx_interactions'] for v in sets.values()]
        title = 'Notx Interactions with Genes by Ratio of KLA {0} to Notx H3K4me2'.format(me2_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by relative amount of H3K4me2', 
                         ylabel='Count of interactions {0}in notx'.format(tss_only and 'with gene TSSs ' or ''),
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        vals = [v['kla_30m_interactions'] for v in sets.values()]
        title = 'KLA 30m Interactions with Genes by Ratio of KLA {0} to Notx H3K4me2'.format(me2_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by relative amount of H3K4me2', 
                         ylabel='Count of interactions {0}in KLA 30m'.format(tss_only and 'with gene TSSs ' or ''),
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        vals = [v['kla_4h_interactions'] for v in sets.values()]
        title = 'KLA 4h Interactions with Genes by Ratio of KLA {0} to Notx H3K4me2'.format(me2_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by relative amount of H3K4me2', 
                         ylabel='Count of interactions {0}in KLA 4h'.format(tss_only and 'with gene TSSs ' or ''),
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)

        vals = [v['kla_30m_over_notx'] for v in sets.values()]
        title = 'Ratio of KLA 30m to Notx Interactions by Ratio of KLA {0} to Notx H3K4me2'.format(me2_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by relative amount of H3K4me2', 
                         ylabel='Ratio of interactions {0}in KLA 30m over notx'.format(tss_only and 'with gene TSSs ' or ''),
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        vals = [v['kla_4h_over_notx'] for v in sets.values()]
        title = 'Ratio of KLA 4h to Notx Interactions by Ratio of KLA {0} to Notx H3K4me2'.format(me2_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by relative amount of H3K4me2', 
                         ylabel='Ratio of interactions {0}in KLA 4h over notx'.format(tss_only and 'with gene TSSs ' or ''),
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        
        
        
