'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'piecharts_by_binding')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_gene.txt'))
    data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
    
    # Make sure we have dimethyl
    data = data[data.filter(like='h3k4me2').max(axis=1) > 10]
    data = data[data['minimal_distance'] >= 1000]
    
    #data = yzer.collapse_strands(data)
    
    print data[(data['h3k4me2_tag_count'] == 0) & (data['gr_dex_tag_count'] > 10)]['ucsc_link_nod']
    print data[(data['h3k4me2_tag_count'] == 0) & (data['p65_kla_tag_count'] > 10)]['ucsc_link_nod']
    
    print sum(data['gr_dex_tag_count'] > 10)
    print sum(data['p65_kla_tag_count'] > 10)
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'feature_vectors.txt'))
    transcripts['id'] = transcripts['glass_transcript_id']
    data = data.merge(transcripts, how='left', on='id', suffixes=['','_trans'])
    
    data = data.fillna(0)
    
    transrepressed = data[(data['kla_1_lfc'] >= 1) & (data['dex_over_kla_1_lfc'] <= -.58)]
    not_trans = data[(data['kla_1_lfc'] < 1) | (data['dex_over_kla_1_lfc'] > -.58)]
    
    supersets = (('Not near transrepressed genes', not_trans),('Near transrepressed genes',transrepressed))
    
    # Plot trans versus not
    yzer.piechart([len(d) for d in zip(*supersets)[1]], zip(*supersets)[0],
                 title='Enhancer-like Subsets by state in KLA+Dex', 
                 save_dir=img_dirpath, show_plot=False)
    
    tfs = [('PU.1','pu_1'),('p65','p65'),('GR','gr')]
    contexts = [('DMSO',''),('Dex','dex'),('KLA','kla'),('KLA+Dex','kla_dex')]
    
    for name, dataset in supersets:
        total_for_set = len(dataset)
        for tf_name, tf in tfs:
            # Get count for enhancer elements with this TF at all
            with_tf = dataset[dataset.filter(like=tf).max(axis=1) > 10]
            without_tf = dataset[dataset.filter(like=tf).max(axis=1) <= 10]
            
            # Plot with TF versus not
            yzer.piechart([ len(without_tf), len(with_tf)], 
                          [s.format(tf_name) for s in ('No {0}', 'Has {0}')],
                         title='Enhancer-like Subsets {0}\nby {1} Occupancy in any state'.format(
                                                                                name.lower(), tf_name), 
                         save_dir=img_dirpath, show_plot=False)
            
            # For those that have the TF, what do they look like in KLA+Dex?
            # First, get default state to compare to...
            for context_name, context in contexts:
                default_context = '{0}_tag_count'.format(tf + (context and ('_' + context) or ''))
                if default_context in with_tf.columns: break
            ratio = 1.5
            none_kla_dex = with_tf[with_tf['{0}_kla_dex_tag_count'.format(tf)] <= 10]
            have_kla_dex = with_tf[with_tf['{0}_kla_dex_tag_count'.format(tf)] > 10]
            loss_kla_dex = have_kla_dex[have_kla_dex['{0}_kla_dex_tag_count'.format(tf)]*ratio < have_kla_dex[default_context]]
            gain_kla_dex = have_kla_dex[have_kla_dex['{0}_kla_dex_tag_count'.format(tf)] > ratio*have_kla_dex[default_context]]
            nc_kla_dex = have_kla_dex[(have_kla_dex['{0}_kla_dex_tag_count'.format(tf)] <= ratio*have_kla_dex[default_context])
                                      & (have_kla_dex['{0}_kla_dex_tag_count'.format(tf)]*ratio >= have_kla_dex[default_context])]
            
            counts = [len(d) for d in (none_kla_dex, loss_kla_dex, nc_kla_dex, gain_kla_dex)]
            labels = [s.format(tf_name) for s in ('No {0} in KLA+Dex', 'Loses {0} in KLA+Dex', 
                                                  'No change in {0}', 'Gains {0} in KLA+Dex')]
            
            yzer.piechart(counts, labels,
                         title='Enhancer-like Subsets {0}\nwith {1} in KLA+Dex verus {2}'.format(
                                                                            name.lower(), tf_name, context_name), 
                         save_dir=img_dirpath, show_plot=False)
            
            
