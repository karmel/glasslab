'''
Created on Feb 12, 2013

@author: karmel

Note that set 1 of Rudensky's Foxp3 chip has 2x as many peaks,
but he seems to use set 2 in the paper.
'''
from __future__ import division
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    grapher = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/TReg_enhancers/2013_04_01'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath, 'motifs')
    graph_dirpath = yzer.get_filename(dirpath, 'piecharts')
    
    min_score = 10
    for antibody in ('me2',):
        data = {}
        celltypes = ['treg','naive','th1','th2']
        for celltype in ('treg','th1'):
            data[celltype] = yzer.import_file(yzer.get_filename(dirpath, 
                                    '{0}_{1}_versus_others_with_foxp3.txt'.format(celltype, 
                                                                     antibody))).fillna(0)
                
            # Filter out promoters
            data[celltype] = data[celltype][data[celltype]['tss_id'] == 0]
        
            # Get venn-diagram "only" sets
            others = celltypes[:]
            others.remove(celltype)
            
            data[celltype + '_only'] = data[celltype][data[celltype]\
                                        [[o + '_id' for o in others]].sum(axis=1) == 0]
            
            # Pairwise
            for other in others:
                other_pair = others[:]
                other_pair.remove(other)
                names = sorted([celltype, other])
                if '_'.join(names) not in data:
                    data['_'.join(names)] = data[celltype][(data[celltype]\
                                            [[o + '_id' for o in other_pair]].sum(axis=1) == 0)\
                                            & (data[celltype][other + '_id'] > 0)]
                    
                    # Triplicate
                    for third in other_pair:
                        fourth = other_pair[:]
                        fourth.remove(third)
                        names = sorted([celltype, other, third])
                        if '_'.join(names) not in data:
                            data['_'.join(names)] = data[celltype][(data[celltype]\
                                                    [[c + '_id' for c in names]].min(axis=1) > 0)\
                                                    & (data[celltype][fourth[0] + '_id'] == 0)]
        
        data['all'] = data['treg'][data['treg']\
                                   [[c + '_id' for c in celltypes]].min(axis=1) > 0]
        
        # Special cases
        # Treg and Th1 shared and not shared, regardless of others
        data['treg_th1_shared'] = data['treg'][data['treg']['th1_id'] > 0]
        data['treg_not_th1'] = data['treg'][data['treg']['th1_id'] == 0]
        data['th1_not_treg'] = data['th1'][data['th1']['treg_id'] == 0]
        
        data['treg_naive_shared'] = data['treg'][data['treg']['naive_id'] > 0]
        data['treg_not_naive'] = data['treg'][data['treg']['naive_id'] == 0]
        #data['naive_not_treg'] = data['naive'][data['naive']['treg_id'] == 0]
        
        for k in sorted(data.keys()):
            subset = data[k]
            total = len(subset)
            acetylated = len(subset[subset['ac_id'] > 0])
            foxp3 = len(subset[subset['foxp3_tag_count'] >= min_score])
            both = len(subset[(subset['ac_id'] > 0) & (subset['foxp3_tag_count'] >= min_score)])
            none = len(subset[(subset['ac_id'] == 0) & (subset['foxp3_tag_count'] < min_score)])
            summary = '''{k}:
            Total enh: {0}
            Acetylated in Treg: {1} ({2}%)
            Foxp3 in Treg: {3} ({4}%)
            Both: {5} ({6}%)
            '''.format(total,
                       acetylated, (acetylated/total)*100,
                       foxp3, (foxp3/total)*100,
                       both, (both/total)*100,
                       k=k)
            print summary
        
            # Draw pie for each group, showing % with foxp3, % with ac, and % with both
            relevant_cells = ', '.join([s.title() for s in k.split('_')])
            counts = [acetylated - both, both, foxp3 - both, none]
            grapher.piechart(counts=counts, 
                             labels=['Has H3K27Ac in Tregs', 'Has Both', 'Has FoxP3 in Tregs', 'Has Neither'], 
                             title='FoxP3 and H3K27Ac at Enhancers\nwith H3K4me2 in {}'.format(relevant_cells), 
                             small_legend=False, 
                             colors=['#FFFB97', '#D5F0CB', '#ABE4FF', 'white'],
                             save_dir=graph_dirpath,
                             show_plot=False)
            
        if False:
            data['with_foxp3'] = data['treg'][data['treg']['foxp3_tag_count'] >= min_score]
            data['without_foxp3'] = data['treg'][data['treg']['foxp3_tag_count'] < min_score]
            
            for k in ('with_foxp3', 'without_foxp3'):
                first_peak = 'treg'
                subset = data[k]
                subset['id'] = subset[first_peak + '_id']
                subset['start'] = subset[first_peak + '_start']
                subset['end'] = subset[first_peak + '_end']
                
                yzer.run_homer(subset, first_peak + '_' + k, motif_dirpath,
                           cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])
    
