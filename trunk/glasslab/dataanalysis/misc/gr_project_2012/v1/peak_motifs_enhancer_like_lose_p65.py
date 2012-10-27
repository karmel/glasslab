'''
Created on Sep 7, 2012

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath,'motifs','from_peaks')
    
    all_data = yzer.import_file(yzer.get_filename(motif_dirpath, 'p65_kla_vectors.txt'))
    size = 100
    if True:
        for ratio in (3,2,1.5):
            enhancers = yzer.import_file(yzer.get_filename(dirpath, 'boxplots_non_refseq_by_p65',
                                                   'enhancer_like_lose_p65_{0}x_change_dsg_only.txt'.format(ratio)))
            enhancers['glass_transcript_id'] = enhancers['id']
            
            # Limit to peaks and touching transcripts, then pull out peaks that intersect our enhancer set
            data = all_data[all_data['touches'] == 't']
            data = data.merge(enhancers, how='right', on='glass_transcript_id', suffixes=['','trans'])
            curr_path = yzer.get_and_create_path(motif_dirpath,  
                                                 'enhancer_like_lose_p65', 'ratio_{0}'.format(ratio))
            # Group them after selecting those that we want
            data = data.groupby(['id','chr_name'],as_index=False).mean()
            
            #bg = yzer.get_filename(motif_dirpath,  
            #        'peak_motifs_by_transcript_lfc', 'p65_kla',
            #        'all','all','all','all_regions_for_homer.txt')
            
            yzer.run_homer(data, 'ratio_{0}'.format(ratio), curr_path, 
                              center=False, reverse=False, preceding=False, size=size,
                              bg=None, cpus=7)
        
