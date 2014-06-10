'''
Created on Feb 12, 2013

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/iTreg_enhancers/2014_04_11/'
    dirpath = yzer.get_path(dirpath)
    
    for ab in ('me2', 'ac'):
        condition = 'treg'
        motif_dirpath = yzer.get_filename(dirpath, 'Motifs')
        
        filename = yzer.get_filename(dirpath, 
                        '{}_{}_enhancers.txt'.format(condition, ab))
        
        data = yzer.import_file(filename)
        data = data.fillna(0)
         
        min_thresh = 20
            
        if False:
            subdata = data[data['tag_count'] > min_thresh]
            subdata = subdata[subdata['tag_count(2)'] > min_thresh]
            subdata = subdata[subdata['tag_count(3)'] <= 0]
            subdata = subdata[subdata['tag_count(4)'] <= 0]
            subdir = 'treg_shared_' + ab
            yzer.run_homer(subdata, 
                    subdir, motif_dirpath,
                    cpus=6, center=True, reverse=False, preceding=False, 
                    size=200, length=[8, 10, 12, 15], mock=True)
            output_file = yzer.get_filename(motif_dirpath, subdir, 
                                            subdir + '_enhancers.txt')
            subdata.to_csv(output_file,
                        header=True, index=False, sep='\t')
            
            subdata = data[data['tag_count'] > min_thresh]
            subdata = subdata[subdata['tag_count(2)'] <= 0]
            subdata = subdata[subdata['tag_count(3)'] <= 0]
            subdata = subdata[subdata['tag_count(4)'] <= 0]
            subdir = 'ntreg_only_' + ab
            yzer.run_homer(subdata, 
                    subdir, motif_dirpath,
                    cpus=6, center=True, reverse=False, preceding=False, 
                    size=200, length=[8, 10, 12, 15], mock=True)
            output_file = yzer.get_filename(motif_dirpath, subdir, 
                                            subdir + '_enhancers.txt')
            subdata.to_csv(output_file,
                        header=True, index=False, sep='\t')
        
        if True:
            subdata = data[data['tag_count'] > min_thresh]
            subdata = subdata[subdata['tag_count(2)'] > min_thresh]
            subdata = subdata[subdata['tag_count(3)'] <= 0]
            subdata = subdata[subdata['tag_count(4)'] <= 0]
            
            # Sort by location and remove all enhancers within x bp of another
            # enhnacer.
            
            min_distance = 5000
            subdata = subdata.sort('start')
            
            prev_idx = None
            omit = []
            for idx, row in subdata.iterrows():
                if prev_idx is not None \
                    and row['start'] - subdata[prev_idx] <= min_distance:
                        omit.append(prev_idx)
                        omit.append(idx)
            omit = set(omit)
            subdata.drop(omit)
            
            subdir = 'treg_shared_isolated_' + ab
            yzer.run_homer(subdata, 
                    subdir, motif_dirpath,
                    cpus=6, center=True, reverse=False, preceding=False, 
                    size=200, length=[8, 10, 12, 15], mock=True)
            output_file = yzer.get_filename(motif_dirpath, subdir, 
                                            subdir + '_enhancers.txt')
            subdata.to_csv(output_file,
                        header=True, index=False, sep='\t')
        