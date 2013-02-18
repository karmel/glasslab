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
    dp = yzer.import_file(yzer.get_filename(dirpath, 'dp_with_{0}.txt'.format(peak))).fillna(0)
    naive = yzer.import_file(yzer.get_filename(dirpath, '{0}_with_dp.txt'.format(peak))).fillna(0)
    
    # Filter out promoters
    dp = dp[dp['tss_id'] == 0]
    naive = naive[naive['tss_id'] == 0]
    
    # Get venn-diagram sets for th1/th2
    dp_only = dp[dp['naive_id'] == 0]
    naive_only = naive[naive['dp_id'] == 0]
    shared = naive[naive['dp_id'] > 0]
    
    datasets = [dp_only, naive_only, shared]
    main_peak = ['dp','naive','naive']
    names = ['dp_only_vs_naive','naive_only_vs_dp','shared_dp_naive',]
    print [len(d) for d in datasets]
    
    
    for i, subset in enumerate(datasets):
        subset['id'] = subset['{0}_id'.format(main_peak[i])]
        subset['start'] = subset['{0}_start'.format(main_peak[i])]
        subset['end'] = subset['{0}_end'.format(main_peak[i])]
        
        yzer.run_homer(subset, names[i], motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

    