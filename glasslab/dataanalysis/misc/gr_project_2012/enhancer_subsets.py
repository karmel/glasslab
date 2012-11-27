'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'piecharts_by_binding','by_genes')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_gene.txt'))
    data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
    
    min_tags = 30
    # Make sure we have dimethyl
    data = data[data.filter(like='h3k4me2').max(axis=1) > min_tags]
    data = data[data['minimal_distance'] >= 1000]
    
    #data = yzer.collapse_strands(data)
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    transcripts['nearest_refseq_transcript_id'] = transcripts['id']
    data = data.merge(transcripts, how='left', on='nearest_refseq_transcript_id', suffixes=['','_trans'])
    
    data = data.fillna(0)
    
    transrepressed = data[(data['kla_1_lfc_trans'] >= 1) & (data['dex_over_kla_1_lfc_trans'] <= -.58)]
    not_trans = data[(data['kla_1_lfc_trans'] < 1) | (data['dex_over_kla_1_lfc_trans'] > -.58)]
    up_in_kla = data[(data['kla_1_lfc_trans'] >= 1) & (data['dex_over_kla_1_lfc_trans'] > -.58)]
    
    supersets = (('All', data),
                 ('Not near transrepressed genes', not_trans), 
                 ('Up in KLA', up_in_kla),
                 ('Near transrepressed genes',transrepressed))
    
    # Plot trans versus not
    yzer.piechart([len(d) for d in zip(*supersets[1:])[1]], zip(*supersets[1:])[0],
                 title='Enhancer-like Subsets by state in KLA+Dex', 
                 save_dir=img_dirpath, show_plot=False)
    
    tfs = [('PU.1','pu_1'),('p65','p65'),('GR','gr')]
    contexts = [('DMSO',''),('Dex','dex'),('KLA','kla'),('KLA+Dex','kla_dex')]
    
    for name, dataset in supersets:
        total_for_set = len(dataset)
        for tf_name, tf in tfs:
            # Get count for enhancer elements with this TF at all
            cols = ['{0}_{1}tag_count'.format(tf, c and (c+'_') or '') for _, c in contexts]
            with_tf = dataset[dataset.filter(items=cols).max(axis=1) > min_tags]
            without_tf = dataset[dataset.filter(items=cols).max(axis=1) <= min_tags]
            
            # Plot with TF versus not
            yzer.piechart([ len(without_tf), len(with_tf)], 
                          [s.format(tf_name) for s in ('No {0}', 'Has {0}')],
                         title='Enhancer-like Subsets {0}\nby {1} Occupancy in any state'.format(
                                                                                name.title(), tf_name), 
                         save_dir=img_dirpath, show_plot=False)
            
            # For those that have the TF, what do they look like in KLA+Dex?
            # First, get default state to compare to...
            default_contexts = []
            for context_name, context in contexts:
                if context == 'kla_dex': continue
                default_context = '{0}_tag_count'.format(tf + (context and ('_' + context) or ''))
                if default_context in with_tf.columns: default_contexts.append((context_name,default_context))
            ratio = 1.5
            have_kla_dex = with_tf[with_tf['{0}_kla_dex_tag_count'.format(tf)] > min_tags]
            
            for context_name, def_context in default_contexts:
                only_in_default = with_tf[(with_tf['{0}_kla_dex_tag_count'.format(tf)] <= min_tags) &
                                       (with_tf[def_context] > min_tags)]
                loss_kla_dex = have_kla_dex[have_kla_dex['{0}_kla_dex_tag_count'.format(tf)]*ratio < have_kla_dex[def_context]]
                gain_kla_dex = have_kla_dex[have_kla_dex['{0}_kla_dex_tag_count'.format(tf)] > ratio*have_kla_dex[def_context]]
                nc_kla_dex = have_kla_dex[(have_kla_dex['{0}_kla_dex_tag_count'.format(tf)] <= ratio*have_kla_dex[def_context])
                                          & (have_kla_dex['{0}_kla_dex_tag_count'.format(tf)]*ratio >= have_kla_dex[def_context])]
                
                counts = [len(d) for d in (only_in_default, loss_kla_dex, nc_kla_dex, gain_kla_dex)]
                labels = [s.format(tf_name) for s in ('{{0}} only in {0}'.format(context_name), 'Loses {0} in KLA+Dex', 
                                                      'No change in {0}', 'Gains {0} in KLA+Dex')]
                
                yzer.piechart(counts, labels,
                             title='Enhancer-like Subsets {0}\nwith {1} in any state: KLA+Dex verus {2}'.format(
                                                                                name.title(), tf_name, context_name), 
                             save_dir=img_dirpath, show_plot=False)
                
            # Of those that have TF in KLA+Dex, what proportion also have other TFs?
            others = tfs[:]
            others.remove((tf_name, tf))
            have_all =  have_kla_dex[have_kla_dex[['{0}_kla_dex_tag_count'.format(other) 
                                                   for _, other in others]].min(axis=1) > min_tags]
            have_none = have_kla_dex[have_kla_dex[['{0}_kla_dex_tag_count'.format(other) 
                                                   for _, other in others]].max(axis=1) <= min_tags]
            
            have_one =  have_kla_dex[have_kla_dex[['{0}_kla_dex_tag_count'.format(other) 
                                                   for _, other in others]].min(axis=1) <= min_tags]
            have_0 = have_one[have_one['{0}_kla_dex_tag_count'.format(others[0][1])] > min_tags]
            have_1 = have_one[have_one['{0}_kla_dex_tag_count'.format(others[1][1])] > min_tags]
            
            counts = [len(d) for d in (have_none, have_0, have_1, have_all)]
            labels = [s.format(others[0][0], others[1][0]) for s in ('No {0} or {1} in KLA+Dex', 
                                                                     'Has {0} in KLA+Dex', 
                                                                     'Has {1} in KLA+Dex', 
                                                                     '{0} and {1} in KLA+Dex')]
            yzer.piechart(counts, labels,
                         title='Enhancer-like Subsets {0}\nwith {1} in KLA+Dex: Other TFs'.format(
                                                                            name.title(), tf_name), 
                         save_dir=img_dirpath, show_plot=False)
            