'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from matplotlib import pyplot
from pandas.core.frame import DataFrame

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'counts_by_gene')
    
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
    # Outer join on transcripts 
    data = transcripts.merge(data, how='left', on='nearest_refseq_transcript_id', suffixes=['_trans',''])
    
    data = data.fillna(0)
    
    transrepressed = data[(data['kla_1_lfc_trans'] >= 1) & (data['dex_over_kla_1_lfc_trans'] <= -.58)]
    not_trans = data[(data['kla_1_lfc_trans'] < 1) | (data['dex_over_kla_1_lfc_trans'] > -.58)]
    up_in_kla = data[(data['kla_1_lfc_trans'] >= 1) & (data['dex_over_kla_1_lfc_trans'] > -.58)]
    
    supersets = (('Up in KLA', up_in_kla),
                 ('Near transrepressed genes',transrepressed),
                 ('All', data),
                 ('Not near transrepressed genes', not_trans), 
                 )
    
    tfs = [('PU.1','pu_1'),('p65','p65'),('GR','gr')]
    contexts = [('DMSO',''),('Dex','dex'),('KLA','kla'),('KLA+Dex','kla_dex')]
    
    labels = ['Total enhancers','With GR in KLA+Dex','Tethered', 'Direct competition, favor to GR', 
              'Indirect competition, favor to GR',
                  'Direct competition, favor to p65', 'Directly co-bound without loss', 
                  'Directly bound novel p65 site', 'Indirectly bound novel p65 site', 
                  'Has GR in Dex, no p65', 'Has GR in KLA+Dex only, no p65']
            
    for name, dataset in supersets:
        # Group by gene id
        grouped = dataset.groupby('nearest_refseq_transcript_id')
        
        # For each gene, we want to count 
        # number of enhancers, and number of enhancers of certain types of interest.
        index, all_stats = [], []
        def count_enhancers(group):
            gr_group = group[group['gr_kla_dex_tag_count'] > min_tags]
            
            stats = [len(group), len(gr_group)]
            
            if gr_group > 0:
                # "Tethered" = No GR in Dex; have GR in KLA+Dex; have p65 in KLA; no loss of p65 in KLA+Dex
                tethered = sum((gr_group['gr_dex_tag_count'] <= min_tags) &
                                   (gr_group['gr_kla_dex_tag_count'] > min_tags) &
                                   (gr_group['p65_kla_tag_count'] > min_tags) &
                                   (gr_group['p65_kla_dex_tag_count']*ratio > gr_group['p65_kla_tag_count'])
                                )
                
                # "Direct competition, favor to GR" = Have GR in Dex; have GR in KLA+Dex; have p65 in KLA; loss of p65 in KLA+Dex
                direct_comp_gr = sum((gr_group['gr_dex_tag_count'] > min_tags) &
                                      (gr_group['gr_kla_dex_tag_count'] > min_tags) &
                                      (gr_group['p65_kla_tag_count'] > min_tags) &
                                      (gr_group['p65_kla_dex_tag_count']*ratio <= gr_group['p65_kla_tag_count'])
                                )
                
                # "Indirect competition" = No GR in Dex; have GR in KLA+Dex; have p65 in KLA; loss of p65 in KLA+Dex
                indirect_comp_gr = sum((gr_group['gr_dex_tag_count'] <= min_tags) &
                                        (gr_group['gr_kla_dex_tag_count'] > min_tags) &
                                        (gr_group['p65_kla_tag_count'] > min_tags) &
                                        (gr_group['p65_kla_dex_tag_count']*ratio <= gr_group['p65_kla_tag_count'])
                                )
                
                # "Direct competition, favor to p65" = Have GR in Dex; loss of GR in KLA+Dex; have p65 in KLA; have p65 in KLA+Dex
                direct_comp_p65 = gr_group[(gr_group['gr_dex_tag_count'] > min_tags) &
                                      (gr_group['gr_kla_dex_tag_count']*ratio <= gr_group['gr_dex_tag_count']) &
                                      (gr_group['p65_kla_tag_count'] > min_tags) &
                                      (gr_group['p65_kla_dex_tag_count'] > min_tags)
                                ]
                
                # "Directly co-bound without loss" = Have GR in Dex; no loss of GR in KLA+Dex; have p65 in KLA; no loss of p65 in KLA+Dex
                cobound = sum((gr_group['gr_dex_tag_count'] > min_tags) &
                                      (gr_group['gr_kla_dex_tag_count']*ratio > gr_group['gr_dex_tag_count']) &
                                      (gr_group['p65_kla_tag_count'] > min_tags) &
                                      (gr_group['p65_kla_dex_tag_count']*ratio > gr_group['p65_kla_tag_count'])
                                )
                
                
                # "Novel p65 sites, direct" = Have GR in Dex; have GR in KLA+Dex; no p65 in KLA; have p65 in KLA+Dex
                direct_novel = sum((gr_group['gr_dex_tag_count'] > min_tags) &
                                       (gr_group['gr_kla_dex_tag_count'] > min_tags) &
                                       (gr_group['p65_kla_tag_count'] <= min_tags) &
                                       (gr_group['p65_kla_dex_tag_count'] > min_tags)
                                )
                
                # "Novel p65 sites, indirect" = No GR in Dex; have GR in KLA+Dex; no p65 in KLA; have p65 in KLA+Dex
                indirect_novel = sum((gr_group['gr_dex_tag_count'] <= min_tags) &
                                         (gr_group['gr_kla_dex_tag_count'] > min_tags) &
                                         (gr_group['p65_kla_tag_count'] <= min_tags) &
                                         (gr_group['p65_kla_dex_tag_count'] > min_tags)
                                )
                
                in_dex_no_p65 = sum((gr_group['gr_dex_tag_count'] > min_tags) &
                                         (gr_group['p65_kla_tag_count'] + gr_group['p65_kla_dex_tag_count'] <= min_tags)
                                 )
                
                kla_only_no_p65 = sum((gr_group['gr_dex_tag_count'] <= min_tags) &
                                         (gr_group['gr_kla_dex_tag_count'] > min_tags) &
                                         (gr_group['p65_kla_tag_count'] + gr_group['p65_kla_dex_tag_count'] <= min_tags)
                                 )
                counts = [tethered, direct_comp_gr, indirect_comp_gr, 
                               direct_comp_p65, cobound,
                               direct_novel, indirect_novel,
                               in_dex_no_p65, kla_only_no_p65]
            else: counts = [0]*9
            
            stats = stats + counts
            all_stats.append(dict(zip(labels, stats))), index.append(group.name)
        
        grouped.apply(count_enhancers)
        # There must be a better way to do this group-apply, but I can't make it turn back into a DF...
        enhancer_counts = DataFrame(all_stats, index=index)
        
        spaced_labels = ['\n'.join(map(' '.join,
                               [l.split()[i:i+2] for i in xrange(0,len(l.split()),2)] )) 
                               for l in labels]
        erna_title = 'Enhancers per Gene by enhancer subtype {0}'.format(name)
        ax = yzer.boxplot([enhancer_counts[col] for col in labels], spaced_labels, 
                         title=erna_title, 
                         xlabel='Subset', 
                         ylabel='Count', 
                         show_outliers=False, show_plot=False, wide=True
                         )
        yzer.ylim(ax, -1, 2)
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, erna_title + '.png'))
        yzer.show_plot()