'''
Created on Feb 12, 2013

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/AND_TCR/Analysis/Chip1/Motifs'
    dirpath = yzer.get_path(dirpath)
    
    for ab in ('ac',):
        for peptide in ('K99A','NoPep','PCC'):
        
            pep_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(peptide, ab))
            filename = yzer.get_filename(pep_dirpath, 
                            '{}_{}_enhancers.txt'.format(peptide, ab))
            
            data = yzer.import_file(filename)
            data = data.fillna(0)


            # me2
            if True:
                yzer.run_homer(data, 
                        'all', pep_dirpath,
                        cpus=6, center=True, reverse=False, preceding=False, 
                        size=200, length=[8, 10, 12, 15], mock=True)
                
                yzer.run_homer(data[data['tag_count'] >= 10], 
                        'tag_thresh_10', pep_dirpath,
                        cpus=6, center=True, reverse=False, preceding=False, 
                        size=200, length=[8, 10, 12, 15], mock=True)
            
