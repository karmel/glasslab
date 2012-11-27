'''
Created on Sep 7, 2012

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer
from glasslab.dataanalysis.misc.gr_project_2012.boxplots_redistribution_pairs import get_high_quality_pairs


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath,'motifs','from_peaks')
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'transcript_vectors.txt'))
    transcripts['glass_transcript_id'] = transcripts['id']
    
    size = 200
    if True:
        all_data = yzer.import_file(yzer.get_filename(dirpath, 'redistribution',
                                               'p65_peaks_bigger_in_kla_dex_with_nearby_bigger_kla_peaks.txt'))
    
        all_data = get_high_quality_pairs(all_data, transcripts)
        
        for super_name, data in (('all', all_data,),
                              #('refseq', all_data[(all_data['score'] > 10) & (all_data['has_refseq'] == 1) 
                              #                & (all_data['touches'] == 't') | (all_data['relationship'] == 'is downstream of')],),
                              #('distal', all_data[(all_data['distal'] == 't')],)
                              ):
            for name, dataset in (('all', data,),
                                  ('with_gr_both', data[(data['tag_count_3'] >= 10) & (data['tag_count_4'] >= 10)]),
                                  ('with_gr_either', data[(data['tag_count_3'] >= 10) | (data['tag_count_4'] >= 10)]),
                                  ('with_gr_kla_dex', data[data['tag_count_3'] >= 10]),
                                  ('with_gr_dex', data[data['tag_count_4'] >= 10]),
                                  #('down_in_dex', data[data['dex_1_lfc'] <= -1]),
                                  #('down_in_kla_dex', data[data['kla_dex_1_lfc'] <= -1]),
                                  #('down_in_kla', data[data['kla_1_lfc'] <= -1]),
                                  #('up_in_dex', data[data['dex_1_lfc'] >= 1]),
                                  #('up_in_kla_dex', data[data['kla_dex_1_lfc'] >= 1]),
                                  #('up_in_kla', data[data['kla_1_lfc'] >= 1]),
                                  #('transrepressed', data[(data['kla_1_lfc'] >= 1) & (data['dex_over_kla_1_lfc'] <= -.58)]),
                                  #('up_in_dex_down_in_kla_dex', data[(data['dex_1_lfc'] >= 1) & (data['kla_dex_1_lfc'] - data['dex_1_lfc'] <= -.58)]),
                                  ):
                curr_path = yzer.get_and_create_path(motif_dirpath,  
                                                     'redistributed_pairs', 'high_quality_pairs_vs_kla_bg',
                                                     super_name, name)
                # Group them after selecting those that we want
                dataset = dataset.groupby(['id','chr_name'],as_index=False).mean()
                
                bg = yzer.get_filename(motif_dirpath,  
                        'peak_motifs_by_transcript_lfc', 'p65_kla',
                        'all','all','all','all_regions_for_homer.txt')
                
                yzer.run_homer(dataset, name, curr_path, 
                                  center=True, reverse=False, preceding=False, size=size,
                                  bg=bg, cpus=7)
            
