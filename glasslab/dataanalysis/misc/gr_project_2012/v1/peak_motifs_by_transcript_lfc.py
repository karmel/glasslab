'''
Created on Sep 7, 2012

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/motifs'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath,'from_peaks')
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    transcripts['glass_transcript_id'] = transcripts['id']
    
    for peak_type in ('gr_dex', 'gr_kla_dex', 'p65_kla_dex','p65_kla'):
        size = 100
        if True:
            all_data = yzer.import_file(yzer.get_filename(motif_dirpath,
                                                   '{0}_vectors.txt'.format(peak_type)))
        
            all_data = all_data.merge(transcripts, how='left', on='glass_transcript_id',suffixes=['','trans'])
            all_data = all_data.fillna(0)
            
            for super_name, data in (#('all', all_data,),
                                  #('refseq', all_data[(all_data['score'] > 10) & (all_data['has_refseq'] == 1) 
                                  #                & (all_data['touches'] == 't') | (all_data['relationship'] == 'is downstream of')],),
                                  #('distal', all_data[(all_data['distal'] == 't')],),
                                  ('enhancer_like_v2', all_data[(all_data['refseq'] == 'f') & (all_data['has_infrastructure'] == 0) 
                                                             & (all_data['length'] < 6000)],),
                                  ):
                for name, dataset in (('all', data,),
                                      #('with_other_kla_dex', data[data['tag_count_3'] >= 10]),
                                      #('no_other_kla_dex', data[data['tag_count_3'] < 10]),
                                      #('with_other_kla_or_dex', data[data['tag_count_4'] >= 10]),
                                      ('no_other_kla_or_dex', data[data['tag_count_4'] < 10]),
                                      #('with_pu_1_kla_dex', data[data['tag_count_5'] >= 10]),
                                      ('no_pu_1_kla_dex', data[data['tag_count_5'] < 10]),
                                      ('gt_partner', data[data['tag_count'] > 1.2*data['tag_count_2']]),
                                      #('lt_partner', data[data['tag_count']*1.2 < data['tag_count_2']]),
                                      ('with_partner', data[data['tag_count_2'] >= 10]),
                                      ('no_partner', data[data['tag_count_2'] < 10]),
                                      #('down_in_dex', data[data['dex_1_lfc'] <= -1]),
                                      #('down_in_kla_dex', data[data['kla_dex_1_lfc'] <= -1]),
                                      #('down_in_kla', data[data['kla_1_lfc'] <= -1]),
                                      #('up_in_dex', data[data['dex_1_lfc'] >= 1]),
                                      #('up_in_kla_dex', data[data['kla_dex_1_lfc'] >= 1]),
                                      #('up_in_kla', data[data['kla_1_lfc'] >= 1]),
                                      #('transrepressed', data[(data['kla_1_lfc'] >= 1) & (data['dex_over_kla_1_lfc'] <= -.58)]),
                                      #('up_in_dex_down_in_kla_dex', data[(data['dex_1_lfc'] >= 1) & (data['kla_dex_1_lfc'] - data['dex_1_lfc'] <= -.58)]),
                                      ):
                    # We have multiple copies of peaks if they align to different transcripts
                    parent_path = yzer.get_and_create_path(motif_dirpath,  
                                                         'peak_motifs_by_transcript_lfc',
                                                         peak_type, super_name)
                    curr_path = yzer.get_and_create_path(parent_path, name)
                    
                    # Group them after selecting those that we want
                    dataset = dataset.groupby(['id','chr_name'],as_index=False).mean()
                    
                    if name != 'all': bg = yzer.get_filename(parent_path, 'all','all','all_regions_for_homer.txt')
                    else: bg = None
                    
                    yzer.run_homer(dataset, name, curr_path, 
                                    center=True, reverse=False, preceding=False, size=size,
                                    cpus=6, bg=bg)
                