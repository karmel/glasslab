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
    img_dirpath = yzer.get_and_create_path(dirpath, 'interactions_by_kla_lfc', tss_only and 'genic' or 'all_interactions')
    
    # File generated in novel_me2_sites
    enhancers = yzer.import_file(yzer.get_filename(data_dirpath,
                    'all_enhancers_with_me2_and_{0}interaction_stats.txt'.format(tss_only and 'tss_' or ''))) 
                    
    for kla_timepoint in ('1h',):
        enhancers['me2_ratio'] = nonzero(enhancers['me2_kla_6h_tag_count_2'])/\
                                    nonzero(enhancers['me2_notx_tag_count_2'])
        
            
        sets = OrderedDict()
        sets['2x GRO in KLA {0}'.format(kla_timepoint)] = enhancers[enhancers[kla_col] > 1]
        sets['No change GRO in KLA {0}'.format(kla_timepoint)] = enhancers[enhancers[kla_col].abs() <= 1]
        sets['1/2 GRO in KLA {0}'.format(kla_timepoint)] = enhancers[enhancers[kla_col] < -1]
        
        labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(sets.keys(), sets.values())]
        
        vals = [v[kla_col] for v in sets.values()]
        title = 'H3K4me2 6h Ratio by KLA {0} Log Fold Change for All Enhancers'.format(kla_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='All distal enhancers by relative amount of H3K4me2', 
                         ylabel='KLA 6h H3Kme2/notx H3Kme2', 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        
        # Now limit to interacting enhancers
        for key, val in sets.iteritems(): sets[key] = val[val['total_interactions'] > 0]
        
        labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(sets.keys(), sets.values())]
        
        vals = [v[kla_col] for v in sets.values()]
        title = 'H3K4me2 6h Ratio by KLA {0} Log Fold Change for Interacting Enhancers'.format(kla_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by KLA LFC', 
                         ylabel='KLA 6h H3Kme2/notx H3Kme2', 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        
        vals = [v['total_interactions'] for v in sets.values()]
        title = 'Total Interactions with Genes by KLA {0} LFC'.format(kla_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by KLA LFC', 
                         ylabel='Count of interactions {0}in all conditions'.format(tss_only and 'with gene TSSs ' or ''), 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        
        
        vals = [v['notx_interactions'] for v in sets.values()]
        title = 'Notx Interactions with Genes by KLA {0} LFC'.format(kla_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by KLA LFC', 
                         ylabel='Count of interactions {0}in notx'.format(tss_only and 'with gene TSSs ' or ''), 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        vals = [v['kla_30m_interactions'] for v in sets.values()]
        title = 'KLA 30m Interactions with Genes by KLA {0} LFC'.format(kla_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by KLA LFC', 
                         ylabel='Count of interactions {0}in KLA 30m'.format(tss_only and 'with gene TSSs ' or ''), 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        vals = [v['kla_4h_interactions'] for v in sets.values()]
        title = 'KLA 4h Interactions with Genes by KLA {0} LFC'.format(kla_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by KLA LFC', 
                         ylabel='Count of interactions {0}in KLA 4h'.format(tss_only and 'with gene TSSs ' or ''), 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)

        vals = [v['kla_30m_over_notx'] for v in sets.values()]
        title = 'Ratio of KLA 30m to Notx Interactions by KLA {0} LFC'.format(kla_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by KLA LFC', 
                         ylabel='Ratio of interactions {0}in KLA 30m over notx'.format(tss_only and 'with gene TSSs ' or ''), 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        vals = [v['kla_4h_over_notx'] for v in sets.values()]
        title = 'Ratio of KLA 4h to Notx Interactions by KLA {0} LFC'.format(kla_timepoint)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Interacting distal enhancers by KLA LFC', 
                         ylabel='Ratio of interactions {0}in KLA 4h over notx'.format(tss_only and 'with gene TSSs ' or ''), 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        
        
        
