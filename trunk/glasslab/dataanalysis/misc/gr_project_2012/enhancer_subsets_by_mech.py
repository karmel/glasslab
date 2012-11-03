'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'piecharts_by_mechanism','tag_count_30')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_gene.txt'))
    data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
    
    min_tags = 30
    ratio = 1.5
    # Make sure we have dimethyl
    data = data[data.filter(like='h3k4me2').max(axis=1) > min_tags]
    data = data[data['minimal_distance'] >= 1000]
    
    #data = yzer.collapse_strands(data)
    
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
        
        # Have GR
        dataset = dataset[dataset.filter(like='gr').max(axis=1) > min_tags]
        total_gr = len(dataset)
        yzer.piechart([total_for_set - total_gr, total_gr], ['No GR', 'Has GR'],
                 title='Enhancer-like Subsets {0}\nby GR Occupancy in any state'.format(name.title()), 
                 save_dir=img_dirpath, show_plot=True)
        
        # "Tethered" = No GR in Dex; have GR in KLA+Dex; have p65 in KLA; no loss of p65 in KLA+Dex
        tethered = dataset[(dataset['gr_dex_tag_count'] <= min_tags) &
                           (dataset['gr_kla_dex_tag_count'] > min_tags) &
                           (dataset['p65_kla_tag_count'] > min_tags) &
                           (dataset['p65_kla_dex_tag_count']*ratio > dataset['p65_kla_tag_count'])
                        ]
        
        # "Direct competition, favor to GR" = Have GR in Dex; have GR in KLA+Dex; have p65 in KLA; loss of p65 in KLA+Dex
        direct_comp_gr = dataset[(dataset['gr_dex_tag_count'] > min_tags) &
                              (dataset['gr_kla_dex_tag_count'] > min_tags) &
                              (dataset['p65_kla_tag_count'] > min_tags) &
                              (dataset['p65_kla_dex_tag_count']*ratio <= dataset['p65_kla_tag_count'])
                        ]
        
        # "Indirect competition" = No GR in Dex; have GR in KLA+Dex; have p65 in KLA; loss of p65 in KLA+Dex
        indirect_comp_gr = dataset[(dataset['gr_dex_tag_count'] <= min_tags) &
                                (dataset['gr_kla_dex_tag_count'] > min_tags) &
                                (dataset['p65_kla_tag_count'] > min_tags) &
                                (dataset['p65_kla_dex_tag_count']*ratio <= dataset['p65_kla_tag_count'])
                        ]
        
        # "Direct competition, favor to p65" = Have GR in Dex; loss of GR in KLA+Dex; have p65 in KLA; have p65 in KLA+Dex
        direct_comp_p65 = dataset[(dataset['gr_dex_tag_count'] > min_tags) &
                              (dataset['gr_kla_dex_tag_count']*ratio <= dataset['gr_dex_tag_count']) &
                              (dataset['p65_kla_tag_count'] > min_tags) &
                              (dataset['p65_kla_dex_tag_count'] > min_tags)
                        ]
        
        # "Directly co-bound without loss" = Have GR in Dex; no loss of GR in KLA+Dex; have p65 in KLA; no loss of p65 in KLA+Dex
        cobound = dataset[(dataset['gr_dex_tag_count'] > min_tags) &
                              (dataset['gr_kla_dex_tag_count']*ratio > dataset['gr_dex_tag_count']) &
                              (dataset['p65_kla_tag_count'] > min_tags) &
                              (dataset['p65_kla_dex_tag_count']*ratio > dataset['p65_kla_tag_count'])
                        ]
        
        
        # "Novel p65 sites, direct" = Have GR in Dex; have GR in KLA+Dex; no p65 in KLA; have p65 in KLA+Dex
        direct_novel = dataset[(dataset['gr_dex_tag_count'] > min_tags) &
                               (dataset['gr_kla_dex_tag_count'] > min_tags) &
                               (dataset['p65_kla_tag_count'] <= min_tags) &
                               (dataset['p65_kla_dex_tag_count'] > min_tags)
                        ]
        
        # "Novel p65 sites, indirect" = No GR in Dex; have GR in KLA+Dex; no p65 in KLA; have p65 in KLA+Dex
        indirect_novel = dataset[(dataset['gr_dex_tag_count'] <= min_tags) &
                                 (dataset['gr_kla_dex_tag_count'] > min_tags) &
                                 (dataset['p65_kla_tag_count'] <= min_tags) &
                                 (dataset['p65_kla_dex_tag_count'] > min_tags)
                        ]
        
        in_dex_no_p65 = dataset[(dataset['gr_dex_tag_count'] > min_tags) &
                                 (dataset.filter(like='p65').max(axis=1) <= min_tags)
                         ]
        
        kla_only_no_p65 = dataset[(dataset['gr_dex_tag_count'] <= min_tags) &
                                 (dataset['gr_kla_dex_tag_count'] > min_tags) &
                                 (dataset.filter(like='p65').max(axis=1) <= min_tags)
                         ]
        
        counts = [len(d) for d in (tethered, direct_comp_gr, indirect_comp_gr, 
                                   direct_comp_p65, cobound,
                                   direct_novel, indirect_novel,
                                   in_dex_no_p65, kla_only_no_p65)]
        counts.append(total_gr - sum(counts))
        
        labels = ['Tethered', 'Direct competition, favor to GR', 'Indirect competition, favor to GR',
                  'Direct competition, favor to p65', 'Directly co-bound without loss', 
                  'Directly bound novel p65 site', 'Indirectly bound novel p65 site', 
                  'Has GR in Dex, no p65', 'Has GR in KLA+Dex only, no p65',
                  'Other']
        yzer.piechart(counts, labels,
                 title='Enhancer-like Subsets {0} with GR\nby Putative Mechanism'.format(name.title()), 
                 small_legend=True,
                 save_dir=img_dirpath, show_plot=True)
        
    