'''
Created on Jun 27, 2014

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/AND_TCR/Analysis/Chips1_2_3/Motifs'
    dirpath = yzer.get_path(dirpath)
    
    for ab in ('me2', 'ac'):
        datasets = []
        for peptide in ('K99A','NoPep','PCC'):
            pep_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(peptide, ab))
            
            
            filename = yzer.get_filename(pep_dirpath, 
                            '{}_{}_enhancers.txt'.format(peptide, ab))
            
            data = yzer.import_file(filename)
            data = data.fillna(0)
            
            # Only in all three reps
            #data = data[(data['p1_id'] > 0) & (data['p2_id'] > 0) & (data['p3_id'] > 0)]
            
            datasets.append(data)
            # ONLY this pep:
            only = sum((data['id(2)'] == 0) & (data['id(3)'] == 0)) 
            print('{peptide} {ab} only: {0}'.format(only, 
                                    peptide=peptide, ab=ab))
        
        # Using PCC, get shared groups:
        data = datasets[2]
        pcc_with_k99a = sum((data['k99a_tag_count'] > 0) 
                            & (data['no_pep_tag_count'] == 0))
        pcc_with_no_pep = sum((data['k99a_tag_count'] == 0) 
                              & (data['no_pep_tag_count'] > 0))
        all_three = sum((data['k99a_tag_count'] > 0) 
                          & (data['no_pep_tag_count'] > 0))
        
        # Finally, get K99A and NoPep
        data = datasets[0]
        k99a_with_no_pep = sum((data['pcc_tag_count'] == 0) 
                               & (data['no_pep_tag_count'] > 0))
        
        print('PCC and K99A {ab}: {0}'.format(pcc_with_k99a, 
                                    peptide=peptide, ab=ab))
        
        print('PCC and NoPep {ab}: {0}'.format(pcc_with_no_pep, 
                                    peptide=peptide, ab=ab))
        
        print('K99A and NoPep {ab}: {0}'.format(k99a_with_no_pep, 
                                    peptide=peptide, ab=ab))
        
        print('All three {ab}: {0}'.format(all_three, 
                                    peptide=peptide, ab=ab))
        