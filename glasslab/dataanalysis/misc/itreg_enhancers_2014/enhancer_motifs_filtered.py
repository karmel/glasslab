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
                
                min_thresh = 20
                data = data[data['tag_count'] > min_thresh]
                yzer.run_homer(data, 
                        'filtered_{}'.format(min_thresh), cond_dirpath,
                        cpus=6, center=True, reverse=False, preceding=False, 
                        size=200, length=[8, 10, 12, 15], mock=True)
            
            if True:
                filename = yzer.get_filename(cond_dirpath, 
                                '{}_{}_enhancers.txt'.format(condition, ab))
                
                data = yzer.import_file(filename)
                data = data.fillna(0)
                
                min_thresh = 20
                subdata = data[data['tag_count'] > min_thresh]
                subdata = subdata[subdata['tag_count(2)'] <= min_thresh]
                subdata = subdata[subdata['tag_count(3)'] <= min_thresh]
                yzer.run_homer(subdata, 
                        'only_filtered_{}'.format(min_thresh), cond_dirpath,
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
                
                min_thresh = 20
                data = data[data['tag_count'] > min_thresh]
                
                both = data[data['tag_count(2)'] > min_thresh]
                only = data[data['tag_count(2)'] <= min_thresh]
                if True:
                    if condition == 'itreg':
                        yzer.run_homer(both, 
                                'both_tregs_filtered_{}'.format(min_thresh), cond_dirpath,
                                cpus=6, center=True, reverse=False, preceding=False, 
                                size=200, length=[8, 10, 12, 15], mock=True)
                    yzer.run_homer(only, 
                            'not_other_treg_filtered_{}'.format(min_thresh), cond_dirpath,
                            cpus=6, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12, 15], mock=True)
                    
        # Motifs for Treg shared by Activated, not iTreg (== CD28?)
        if False:
            condition = 'treg'
            filename = yzer.get_filename(cond_dirpath, 
                                '{}_{}_enhancers.txt'.format(condition, ab))
                
            data = yzer.import_file(filename)
            data = data.fillna(0)
            
            min_thresh = 20
            data = data[data['tag_count'] > min_thresh]
            
            with_act = data[(data['tag_count(2)'] <= min_thresh) &
                            (data['tag_count(3)'] > min_thresh)]
            yzer.run_homer(with_act, 
                    'with_act_not_itreg_filtered_{}'.format(min_thresh), cond_dirpath,
                    cpus=6, center=True, reverse=False, preceding=False, 
                    size=200, length=[8, 10, 12, 15], mock=True)
        
        # Motifs for Treg shared by iTreg, not Activated (== Foxp3-driven?)
        if False:
            condition = 'treg'
            filename = yzer.get_filename(cond_dirpath, 
                                '{}_{}_enhancers.txt'.format(condition, ab))
                
            data = yzer.import_file(filename)
            data = data.fillna(0)
            
            min_thresh = 20
            data = data[data['tag_count'] > min_thresh]
            
            with_itreg = data[(data['tag_count(2)'] > min_thresh) &
                              (data['tag_count(3)'] <= min_thresh)]
            
            yzer.run_homer(with_itreg, 
                    'with_itreg_not_act_filtered_{}'.format(min_thresh), cond_dirpath,
                    cpus=6, center=True, reverse=False, preceding=False, 
                    size=200, length=[8, 10, 12, 15], mock=True)
        
        # Motifs for all three shared
        if False:
            condition = 'treg'
            filename = yzer.get_filename(cond_dirpath, 
                                '{}_{}_enhancers.txt'.format(condition, ab))
                
            data = yzer.import_file(filename)
            data = data.fillna(0)
            
            min_thresh = 20
            data = data[data['tag_count'] > min_thresh]
            
            with_both = data[(data['tag_count(3)'] > min_thresh) &
                            (data['tag_count(2)'] > min_thresh)]
            yzer.run_homer(with_both, 
                    'all_three_filtered_{}'.format(min_thresh), cond_dirpath,
                    cpus=6, center=True, reverse=False, preceding=False, 
                    size=200, length=[8, 10, 12, 15], mock=True)
                    