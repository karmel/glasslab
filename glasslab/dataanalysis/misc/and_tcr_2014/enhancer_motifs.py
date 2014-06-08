'''
Created on Feb 12, 2013

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/AND_TCR/Analysis/Chips1_2/Motifs'
    dirpath = yzer.get_path(dirpath)
    
    for ab in ('me2', 'ac'):
        for peptide in ('K99A','NoPep','PCC'):
            pep_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(peptide, ab))
            
            if False:
                filename = yzer.get_filename(pep_dirpath, 
                                '{}_{}_enhancers.txt'.format(peptide, ab))
                
                data = yzer.import_file(filename)
                data = data.fillna(0)
                
    
                yzer.run_homer(data, 
                        'all', pep_dirpath,
                        cpus=6, center=True, reverse=False, preceding=False, 
                        size=200, length=[8, 10, 12, 15], mock=True)
                
                yzer.run_homer(data[data['tag_count'] >= 10], 
                        'tag_thresh_10', pep_dirpath,
                        cpus=6, center=True, reverse=False, preceding=False, 
                        size=200, length=[8, 10, 12, 15], mock=True)
                
            if False:
                filename = yzer.get_filename(pep_dirpath, 
                                '{}_{}_only_enhancers.txt'.format(peptide, ab))
                
                data = yzer.import_file(filename)
                data = data.fillna(0)
    
                if True:
                    yzer.run_homer(data, 
                            'only', pep_dirpath,
                            cpus=6, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12, 15], mock=True)
                    
            if False and peptide == 'PCC' and ab == 'me2':
                filename = yzer.get_filename(pep_dirpath, 
                                '{}_{}_enhancers_batf.txt'.format(peptide, ab))
                
                data = yzer.import_file(filename)
                data = data.fillna(0)
    

                if True:
                    subset = data[data['batf_tag_count'] > 0]
                    yzer.run_homer(subset, 
                            'with_batf', pep_dirpath,
                            cpus=6, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12, 15], mock=True)
                    
                    subset = data[data['irf4_tag_count'] > 0]
                    yzer.run_homer(subset, 
                            'with_irf4', pep_dirpath,
                            cpus=6, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12, 15], mock=True)
                    
                    subset = subset[subset['batf_tag_count'] > 0]
                    yzer.run_homer(subset, 
                            'with_batf_irf4', pep_dirpath,
                            cpus=6, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12, 15], mock=True)

            if False and peptide == 'PCC':
                # de novo shared by PCC and K99A
                filename = yzer.get_filename(pep_dirpath, 
                                '{}_{}_enhancers_batf.txt'.format(peptide, ab))
                
                data = yzer.import_file(filename)
                data = data.fillna(0)
    
                subset = data[data['k99a_tag_count'] > 0]
                subset = subset[subset['no_pep_tag_count'] == 0]
                
                yzer.run_homer(subset, 
                            'with_k99a_only', pep_dirpath,
                            cpus=6, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12, 15], mock=True)
                    
            
            if True and peptide != 'NoPep':
                # Newly active PCC and K99A
                filename = yzer.get_filename(pep_dirpath, 
                                '{}_{}_enhancers_batf.txt'.format(peptide, ab))
                
                data = yzer.import_file(filename)
                data = data.fillna(0)
    
                subset = data[data['no_pep_tag_count'] == 0]
                
                yzer.run_homer(subset, 
                            'not_no_pep', pep_dirpath,
                            cpus=6, center=True, reverse=False, preceding=False, 
                            size=200, length=[8, 10, 12, 15], mock=True)
                    