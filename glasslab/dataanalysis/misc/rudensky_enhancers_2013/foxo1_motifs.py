'''
Created on Feb 20, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer

if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/Foxo1'
    dirpath = yzer.get_path(dirpath)
    motifs_dirpath = yzer.get_and_create_path(dirpath, 'motifs')
    
    peak_pretty = 'Foxo1'
    peak = peak_pretty.lower()
    foxo1 = yzer.import_file(yzer.get_filename(dirpath, '{0}_with_foxp3.txt'.format(peak))).fillna(0)
    
    datasets = [('foxo1_all', foxo1),
                ('foxo1_tss', foxo1[foxo1['tss_id'] > 0]),
                ('foxo1_enhancers', foxo1[foxo1['tss_id'] == 0]),
                ]
    for name, subset in datasets:
        subset['id'] = subset['{0}_id'.format(peak)]
        subset['start'] = subset['{0}_start'.format(peak)]
        subset['end'] = subset['{0}_end'.format(peak)]
        
        yzer.run_homer(subset, name, motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

    