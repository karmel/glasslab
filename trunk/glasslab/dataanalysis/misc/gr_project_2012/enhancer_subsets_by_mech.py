'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from matplotlib import pyplot

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'piecharts_by_mechanism','by_genes')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_gene.txt'))
    data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
    
    draw_pies = True
    min_tags = 30
    ratio = 1.5
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
    if draw_pies:
        yzer.piechart([len(d) for d in zip(*supersets[1:])[1]], zip(*supersets[1:])[0],
                     title='Enhancer-like Subsets by state in KLA+Dex', 
                     save_dir=img_dirpath, show_plot=False)
        
    tfs = [('PU.1','pu_1'),('p65','p65'),('GR','gr')]
    contexts = [('DMSO',''),('Dex','dex'),('KLA','kla'),('KLA+Dex','kla_dex')]
    
    for name, dataset in supersets:
        total_for_set = len(dataset)
        
        # Have GR
        dataset = dataset[dataset['gr_dex_tag_count'] + dataset['gr_kla_dex_tag_count'] > min_tags]
        total_gr = len(dataset)
        if draw_pies:
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
                                 (dataset['p65_kla_tag_count'] + dataset['p65_kla_dex_tag_count'] <= min_tags)
                         ]
        
        kla_only_no_p65 = dataset[(dataset['gr_dex_tag_count'] <= min_tags) &
                                 (dataset['gr_kla_dex_tag_count'] > min_tags) &
                                 (dataset['p65_kla_tag_count'] + dataset['p65_kla_dex_tag_count'] <= min_tags)
                         ]
        
        sets = [tethered, direct_comp_gr, indirect_comp_gr, 
                           direct_comp_p65, cobound,
                           direct_novel, indirect_novel,
                           in_dex_no_p65, kla_only_no_p65]
        counts = [len(d) for d in sets]
        counts.append(total_gr - sum(counts))
        
        labels = ['Tethered', 'Direct competition, favor to GR', 'Indirect competition, favor to GR',
                  'Direct competition, favor to p65', 'Directly co-bound without loss', 
                  'Directly bound novel p65 site', 'Indirectly bound novel p65 site', 
                  'Has GR in Dex, no p65', 'Has GR in KLA+Dex only, no p65',
                  'Other']
        if draw_pies:
            yzer.piechart(counts, labels,
                     title='Enhancer-like Subsets {0} with GR\nby Putative Mechanism'.format(name.title()), 
                     small_legend=True,
                     save_dir=img_dirpath, show_plot=True)
        
        enh_change = [d['dex_over_kla_1_lfc'] for d in [data, dataset] + sets]
        change_labels = ['All Enhancer-like', 'All {0}'.format(name)] + labels
        change_labels = ['\n'.join(map(' '.join,
                               [l.split()[i:i+2] for i in xrange(0,len(l.split()),2)] )) 
                               for l in change_labels]
        erna_title = 'eRNA Change by enhancer subtype {0}- KLA+Dex versus KLA'.format(name)
        ax = yzer.boxplot(enh_change, change_labels, 
                         title=erna_title, 
                         xlabel='Subset', 
                         ylabel='log2(KLA+Dex GRO-seq/KLA GRO-seq)', 
                         show_outliers=False, show_plot=False, wide=True
                         )
        yzer.ylim(ax, -4, 2)
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, erna_title + '.png'))
        yzer.show_plot()