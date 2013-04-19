'''
Created on Apr 18, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/TReg_enhancers/2013_04_01'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath, 'motifs')

    foxp3 = yzer.import_file(yzer.get_filename(dirpath, 
                            'foxp3_with_stat1.txt')).fillna(0)
    stat1 = yzer.import_file(yzer.get_filename(dirpath, 
                            'stat1_with_foxp3.txt')).fillna(0)
                            
    print len(foxp3), len(stat1)
    print sum(foxp3['stat1_id'] > 0)
    print sum(foxp3['stat1_id'] > 0)/len(foxp3)
    print sum(stat1['foxp3_id'] > 0)
    print sum(stat1['foxp3_id'] > 0)/len(stat1)
    
    foxp3_enh = foxp3[(foxp3['tss_me2_id'] == 0) & (foxp3['tss_id'] == 0)]
    foxp3_tss = foxp3[(foxp3['tss_me2_id'] > 0) | (foxp3['tss_id'] > 0)]
    print len(foxp3_enh)
    print sum(foxp3_enh['stat1_id'] > 0)/len(foxp3_enh)
    print len(foxp3_tss)
    print sum(foxp3_tss['stat1_id'] > 0)/len(foxp3_tss)
    
    foxp3_with_stat = foxp3[foxp3['stat1_id'] > 0]
    if False:
        grapher = SeqGrapher()
        grapher.scatterplot(foxp3_with_stat,
                            xcolname='foxp3_tag_count',
                            ycolname='stat1_tag_count',
                            log=True, 
                            show_plot=True)
        
    if False:
        subsets = [('all', foxp3_with_stat),
                   ('enh', foxp3_enh[foxp3_enh['stat1_id'] > 0]),
                   ('tss', foxp3_tss[foxp3_tss['stat1_id'] > 0]),
                   ]
        for k, subset in subsets:
            first_peak = 'foxp3'
            subset['id'] = subset[first_peak + '_id']
            subset['start'] = subset[first_peak + '_start']
            subset['end'] = subset[first_peak + '_end']
            
            yzer.run_homer(subset, first_peak + '_with_stat1_' + k, motif_dirpath,
                           cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

    if True:
        subsets = [('all', foxp3[foxp3['stat1_id'] == 0]),
                   ('enh', foxp3_enh[foxp3_enh['stat1_id'] == 0]),
                   ('tss', foxp3_tss[foxp3_tss['stat1_id'] == 0]),
                   ]
        for k, subset in subsets:
            first_peak = 'foxp3'
            subset['id'] = subset[first_peak + '_id']
            subset['start'] = subset[first_peak + '_start']
            subset['end'] = subset[first_peak + '_end']
            
            yzer.run_homer(subset, first_peak + '_without_stat1_' + k, motif_dirpath,
                           cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

