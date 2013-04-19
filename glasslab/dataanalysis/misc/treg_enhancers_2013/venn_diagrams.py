'''
Created on Feb 12, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/TReg_enhancers/2013_04_01'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath, 'motifs')
    
    
    for antibody in ('me2',):
        data = {}
        celltypes = ['treg','naive','th1','th2']
        for celltype in celltypes:
            data[celltype] = yzer.import_file(yzer.get_filename(dirpath, 
                                    '{0}_{1}_versus_others.txt'.format(celltype, 
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
        
        for k in sorted(data.keys()):
            subset = data[k] 
            print k, len(subset)
        
        if True:
            for k in sorted(data.keys(), reverse=True):
                subset = data[k]    
                if k in celltypes or len(subset) < 1000: continue
                
                if k in ('treg_th1_shared','treg_not_th1','th1_not_treg'): continue
                first_peak = k=='all' and 'naive' or k.split('_')[0]
                subset['id'] = subset[first_peak + '_id']
                subset['start'] = subset[first_peak + '_start']
                subset['end'] = subset[first_peak + '_end']
                
                yzer.run_homer(subset, 'four_way_venn_' + k, motif_dirpath,
                           cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])
    
