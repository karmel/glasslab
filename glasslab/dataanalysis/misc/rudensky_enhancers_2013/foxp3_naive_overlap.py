'''
Created on Feb 20, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer

if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/Rudensky_enhancers'
    dirpath = yzer.get_path(dirpath)
    motifs_dirpath = yzer.get_and_create_path(dirpath, 'motifs')
    
    peak_pretty = 'Foxp3'
    peak = peak_pretty.lower()
    foxp3 = yzer.import_file(yzer.get_filename(dirpath, '{0}_1_with_naive_me2.txt'.format(peak))).fillna(0)
    naive = yzer.import_file(yzer.get_filename(dirpath, 'naive_me2_with_{0}.txt'.format(peak))).fillna(0)
    
    # Filter out promoters
    foxp3 = foxp3[foxp3['tss_id'] == 0]
    naive = naive[naive['tss_id'] == 0]
    
    # Get venn-diagram sets for foxp3/me2
    only_foxp3 = foxp3[foxp3['naive_id'] == 0]
    only_naive = naive[naive['foxp3_1_id'] == 0]
    shared = foxp3[foxp3['naive_id'] > 0]
    print len(only_foxp3), len(only_naive), len(shared)
    
    # Now factor in foxo1
    foxp3_no_foxo1 = foxp3[foxp3['foxo1_id'] == 0]
    foxp3_with_foxo1 = foxp3[foxp3['foxo1_id'] > 0]
    only_foxp3_no_foxo1 = only_foxp3[only_foxp3['foxo1_id'] == 0]
    only_foxp3_with_foxo1 = only_foxp3[only_foxp3['foxo1_id'] > 0]
    shared_no_foxo1 = shared[shared['foxo1_id'] == 0]
    shared_with_foxo1 = shared[shared['foxo1_id'] > 0]
    
    datasets = [foxp3, only_foxp3_no_foxo1, only_foxp3_with_foxo1, shared_no_foxo1,shared_with_foxo1,
                foxp3_no_foxo1, foxp3_with_foxo1]
    main_peak = ['foxp3_1']*len(datasets)
    names = ['all_foxp3','only_foxp3_no_foxo1','only_foxp3_with_foxo1','shared_no_foxo1','shared_with_foxo1',
             'foxp3_no_foxo1','foxp3_with_foxo1']
    print [len(d) for d in datasets]
    
    
    for i, subset in enumerate(datasets):
        if i < 5: continue
        subset['id'] = subset['{0}_id'.format(main_peak[i])]
        subset['start'] = subset['{0}_start'.format(main_peak[i])]
        subset['end'] = subset['{0}_end'.format(main_peak[i])]
        
        yzer.run_homer(subset, names[i], motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

    