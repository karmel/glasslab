'''
Created on Mar 1, 2014

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/iTreg_enhancers/2014_02_14/Motifs'
    dirpath = yzer.get_path(dirpath)
    
    for ab in ('me2', 'ac'):
        for condition in ('treg','itreg'):
            cond_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(condition, ab))
    
            filename = yzer.get_filename(cond_dirpath, 
                    '{}_{}_enhancers.txt'.format(condition, ab))
                    
            data = yzer.import_file(filename)
            data = data.fillna(0)
            
            print('Total {} {}'.format(condition,ab))
            print(len(data))
    
            min_thresh = 20
            if True:
                data = data[data['tag_count'] > min_thresh]
                print('Filtered {} {}'.format(condition,ab))
                print(len(data))
            
            subdata = data[(data['tag_count(2)'] > min_thresh) &
                           (data['tag_count(3)'] > min_thresh)
                           ]
            print('All three ' + ab)
            print(len(subdata), len(subdata)/len(data))
            
            subdata = data[(data['tag_count(2)'] > min_thresh) &
                           (data['tag_count(3)'] <= min_thresh)
                           ]
            print('{} and p2 {}'.format(condition,ab))
            print(len(subdata), len(subdata)/len(data))
            
            subdata = data[(data['tag_count(2)'] <= min_thresh) &
                           (data['tag_count(3)'] > min_thresh)
                           ]
            print('{} and p3 {}'.format(condition,ab))
            print(len(subdata), len(subdata)/len(data))
            
            
            subdata = data[(data['tag_count(2)'] <= min_thresh) &
                           (data['tag_count(3)'] <= min_thresh)
                           ]
            print('{} only {}'.format(condition,ab))
            print(len(subdata), len(subdata)/len(data))
            
            
            