'''
Created on Feb 12, 2013

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/iTreg_enhancers/2014_02_14/Motifs'
    dirpath = yzer.get_path(dirpath)
    
    for ab in ('me2', 'ac'):
        for condition in ('treg','itreg','activated'):
            cond_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(condition, ab))
            
            if False:
                filename = yzer.get_filename(cond_dirpath, 
                                '{}_{}_enhancers.txt'.format(condition, ab))
                
                data = yzer.import_file(filename)
                data = data.fillna(0)
                
                min_thresh = 0
                cutoff = 10000
                data = data[data['tag_count'] > min_thresh]
                data = data.sort('tag_count', ascending=False)[:cutoff]
                
                yzer.run_homer(data, 
                        'top_{}'.format(cutoff), cond_dirpath,
                        cpus=6, center=True, reverse=False, preceding=False, 
                        size=200, length=[8, 10, 12, 15], mock=True)
        
        for condition in ('itreg','treg'):
            cond_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(condition, ab))
            if False:
                filename = yzer.get_filename(cond_dirpath, 
                                '{}_{}_enhancers.txt'.format(condition, ab))
                
                data = yzer.import_file(filename)
                data = data.fillna(0)
                
                min_thresh = 0
                cutoff = 10000
                data = data[data['tag_count'] > min_thresh]
                data = data.sort('tag_count', ascending=False)[:cutoff]
                
                
                both = data[data['tag_count(2)'] > min_thresh]
                only = data[data['tag_count(2)'] <= min_thresh]
                if True:
                    yzer.run_homer(both, 
                            'both_tregs_top_{}'.format(cutoff), cond_dirpath,
                            cpus=6, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12, 15], mock=True)
                    yzer.run_homer(only, 
                            'not_other_treg_top_{}'.format(cutoff), cond_dirpath,
                            cpus=6, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12, 15], mock=True)
                    
        