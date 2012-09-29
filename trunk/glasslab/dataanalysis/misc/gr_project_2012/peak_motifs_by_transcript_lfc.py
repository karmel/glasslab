'''
Created on Sep 7, 2012

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/motifs/from_peaks/vs_non_dsg'
    dirpath = yzer.get_path(dirpath)
    
    for peak_type in ('gr_dex', 'gr_kla_dex','gr_dex_dsg', 'gr_kla_dex_dsg',):
        size = 200
        if True:
            all_data = yzer.import_file(yzer.get_filename(dirpath,
                                                   '{0}_vectors.txt'.format(peak_type)))
        
            all_data = all_data.fillna(0)
            
            for super_name, data in (('all', all_data,),
                                  ('refseq', all_data[(all_data['score'] > 10) & (all_data['has_refseq'] == 1) 
                                                  & (all_data['touches'] == 't') | (all_data['relationship'] == 'is downstream of')],),
                                  ('distal', all_data[(all_data['distal'] == 't')],)):
                for name, dataset in (('all', data,),
                                      ):
                    # We have multiple copies of peaks if they align to different transcripts
                    curr_path = yzer.get_and_create_path(dirpath,  
                                                         peak_type, super_name, name)
                    # Group them after selecting those that we want
                    dataset = dataset.groupby(['id','chr_name'],as_index=False).mean()
                    yzer.run_homer(None, name, curr_path, 
                                  center=True, reverse=False, preceding=False, size=size,
                                  cpus=6, bg=bg,
                                  files_already_prepped=True)
                
