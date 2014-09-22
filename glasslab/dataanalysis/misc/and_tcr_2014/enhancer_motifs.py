'''
Created on Feb 12, 2013

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/AND_TCR/Analysis/IDR/Motifs'
    dirpath = yzer.get_path(dirpath)
    
    for ab in ('me2', ): #'ac'):
        for peptide in ('K99A','NoPep','PCC'):
            pep_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(peptide, ab))
            
            if True:
                filename = yzer.get_filename(pep_dirpath, 
                                '{}_{}_enhancers.txt'.format(peptide, ab))
                
                data = yzer.import_file(filename)
                data = data.fillna(0)
                
    
                if True:
                    yzer.run_homer(data, 
                        'all', pep_dirpath,
                        cpus=8, center=True, reverse=False, preceding=False, 
                        size=200, length=[8, 10, 12], mock=True)
                
        
                if True:
                    subset = data[(data['id(2)'] == 0) & (data['id(3)'] == 0)]
                    yzer.run_homer(subset, 
                            'only', pep_dirpath,
                            cpus=8, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12], mock=True)
                    
                if False and peptide == 'PCC' and ab == 'me2':
                    subset = data[data['batf_tag_count'] > 0]
                    yzer.run_homer(subset, 
                            'with_batf', pep_dirpath,
                            cpus=8, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12], mock=True)
                    
                    subset = data[data['irf4_tag_count'] > 0]
                    yzer.run_homer(subset, 
                            'with_irf4', pep_dirpath,
                            cpus=8, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12], mock=True)
                    
                    subset = data[(data['irf4_tag_count'] > 0) & (data['batf_tag_count'] > 0)]
                    yzer.run_homer(subset, 
                            'with_batf_irf4', pep_dirpath,
                            cpus=8, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12], mock=True)

                if True and peptide != 'NoPep':
                    # Newly active PCC and K99A
                    subset = data[data['no_pep_tag_count'] == 0]
                    
                    yzer.run_homer(subset, 
                                'not_no_pep', pep_dirpath,
                                cpus=8, center=True, reverse=False, preceding=False, 
                                size=200, length=[8, 10, 12], mock=True)
                    
                if True:
                    # Shared for all
                    subset = data[(data['id(2)'] > 0) & (data['id(3)'] > 0)]
                    yzer.run_homer(subset, 
                            'shared_by_all', pep_dirpath,
                            cpus=8, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12], mock=True)
                        