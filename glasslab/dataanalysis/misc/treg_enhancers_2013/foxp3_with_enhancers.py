'''
Created on Feb 12, 2013

@author: karmel

Note that set 1 of Rudensky's Foxp3 chip has 2x as many peaks,
but he seems to use set 2 in the paper.
'''
from __future__ import division
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer
from collections import OrderedDict


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/TReg_enhancers/2013_04_01'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath, 'motifs')
    
    min_score = 10
    for antibody in ('me2',):
        celltypes = ['treg','naive','th1','th2']
        
        stats = OrderedDict()
        all_data = yzer.import_file(yzer.get_filename(dirpath, 
                                'foxp3_with_treg_enhancers.txt')).fillna(0)
        stats['all'] = len(all_data)
        # Filter out promoters
        data = all_data[all_data['tss_id'] == 0]
        stats['enhancers'] = len(data)
        stats['me2_treg'] = sum(data['me2_id'] > 0)
        stats['ac_treg'] = sum(data['ac_id'] > 0)
        
        stats['me2_treg_naive'] = sum((data['me2_id'] > 0) & (data['naive_id'] > 0))
        stats['me2_treg_th1'] = sum((data['me2_id'] > 0) & (data['th1_id'] > 0))
        stats['me2_treg_naive_th1'] = sum((data['me2_id'] > 0) 
                                          & (data['naive_id'] > 0)
                                          & (data['th1_id'] > 0))
        stats['me2_treg_only'] = sum((data['me2_id'] > 0) 
                                          & (data['naive_id'] == 0)
                                          & (data['th1_id'] == 0))
        stats['me2_naive_only'] = sum((data['me2_id'] == 0) 
                                          & (data['naive_id'] > 0)
                                          & (data['th1_id'] == 0))
        stats['me2_th1_only'] = sum((data['me2_id'] == 0) 
                                          & (data['naive_id'] == 0)
                                          & (data['th1_id'] > 0))
        
        
        print 'Set\tCount\t% of Total\t% of Enh'
        for k, count in stats.iteritems():
            print '{}\t{}\t{:.2f}\t{:.2f}'.format(k,
                                           count,
                                           count/stats['all']*100,
                                           count/stats['enhancers']*100,
                                           )
        
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
    
