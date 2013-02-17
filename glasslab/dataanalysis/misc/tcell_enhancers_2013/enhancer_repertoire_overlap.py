'''
Created on Feb 17, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer

if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/H3K4me2/Analysis/venn_diagrams'
    dirpath = yzer.get_path(dirpath)
    motifs_dirpath = yzer.get_and_create_path(dirpath, 'motifs')
    
    peak_pretty = 'H3K4me2'
    peak = peak_pretty.lower()
    th1 = yzer.import_file(yzer.get_filename(dirpath, 'th1_with_th2_me1_and_{0}.txt'.format(peak))).fillna(0)
    th2 = yzer.import_file(yzer.get_filename(dirpath, 'th2_with_th1_me1_and_{0}.txt'.format(peak))).fillna(0)
    naive = yzer.import_file(yzer.get_filename(dirpath, '{0}_with_th1_and_th2_me1.txt'.format(peak))).fillna(0)
    
    # Filter out promoters
    th1 = th1[th1['tss_id'] == 0]
    th2 = th2[th2['tss_id'] == 0]
    naive = naive[naive['tss_id'] == 0]
    
    # Get venn-diagram sets for th1/th2
    only_th1 = th1[th1['th2_id'] == 0]
    only_th2 = th2[th2['th1_id'] == 0]
    shared = th1[th1['th2_id'] > 0]
    shared_check = th2[th2['th2_id'] > 0]
    print len(only_th1), len(only_th2), len(shared), len(shared_check)
    
    # Now factor in me2
    only_th1_no_me2 = only_th1[only_th1['naive_id'] == 0]
    only_th1_with_me2 = only_th1[only_th1['naive_id'] > 0]
    only_th2_no_me2 = only_th2[only_th2['naive_id'] == 0]
    only_th2_with_me2 = only_th2[only_th2['naive_id'] > 0]
    shared_no_me2 = shared[shared['naive_id'] == 0]
    shared_with_me2 = shared[shared['naive_id'] > 0]
    naive_only = naive[(naive['th1_id'] == 0) & (naive['th2_id'] == 0)]
    
    datasets = [only_th1_no_me2,only_th1_with_me2,only_th2_no_me2,only_th2_with_me2,
                shared_no_me2,shared_with_me2, naive_only]
    main_peak = ['th1','th1','th2','th2','naive','naive','naive']
    names = ['only_th1_no_me2','only_th1_with_me2','only_th2_no_me2','only_th2_with_me2',
                'shared_no_me2','shared_with_me2', 'naive_only']
    print [len(d) for d in datasets]
    
    
    for i, subset in enumerate(datasets):
        subset['id'] = subset['{0}_id'.format(main_peak[i])]
        subset['start'] = subset['{0}_start'.format(main_peak[i])]
        subset['end'] = subset['{0}_end'.format(main_peak[i])]
        
        yzer.run_homer(subset, names[i], motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

    