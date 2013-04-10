'''
Created on Mar 28, 2013

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/TReg_enhancers/2013_04_01'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath, 'motifs')
    
    for antibody in ('me2','ac'):
        all_treg = yzer.import_file(yzer.get_filename(dirpath, 
                                    'treg_with_naive_{0}.txt'.format(antibody))).fillna(0)
        all_naive = yzer.import_file(yzer.get_filename(dirpath,
                                    'naive_with_treg_{0}.txt'.format(antibody))).fillna(0)
        
        # Filter out promoters
        treg = all_treg[all_treg['tss_id'] == 0]
        naive = all_naive[all_naive['tss_id'] == 0]
        
        # Get venn-diagram sets for foxp3/me2
        only_treg = treg[treg['naive_id'] == 0]
        only_naive = naive[naive['treg_id'] == 0]
        shared = treg[treg['naive_id'] > 0]
        print len(only_treg), len(only_naive), len(shared)
        
        datasets = [all_treg, all_naive, treg, naive, only_treg, only_naive, shared]
        main_peak = ['treg','naive','treg','naive','treg','naive','treg']
        names = [x.format(antibody) for x in ('all_treg_{0}_peaks',
                                              'all_naive_{0}_peaks',
                                              'all_treg_{0}_enhancers',
                                              'all_naive_{0}_enhancers',
                                              'only_treg_{0}_enhancers',
                                              'only_naive_{0}_enhancers',
                                              'shared_{0}_enhancers')]
        for i, subset in enumerate(datasets):
            subset['id'] = subset['{0}_id'.format(main_peak[i])]
            subset['start'] = subset['{0}_start'.format(main_peak[i])]
            subset['end'] = subset['{0}_end'.format(main_peak[i])]
            
            yzer.run_homer(subset, names[i], motif_dirpath,
                       cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

    