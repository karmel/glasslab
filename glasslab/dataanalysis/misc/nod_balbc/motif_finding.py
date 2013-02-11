'''
Created on Feb 8, 2013

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    base_dirpath = yzer.get_path('karmel/GlassLab/Notes_and_Reports/NOD_BALBc/ThioMacs/Analysis_2013_02/')
    dirpath = yzer.get_and_create_path(base_dirpath, 'motifs/')
    filename = yzer.get_filename(base_dirpath, 'transcript_vectors.txt')
    data = yzer.import_file(filename)
    data = data.fillna(0)
    
    # Promoters
    if False:
        refseq = data[data['has_refseq'] == 1]
        refseq = refseq[refseq['transcript_score'] >= 4]
        if True:
            yzer.run_homer(refseq, 'refseq_promoter', dirpath,
                       cpus=6, center=False, reverse=False, preceding=True, size=400, length=[8, 10, 12, 15])

        bg = yzer.get_filename(dirpath, 'refseq_promoter/refseq_promoter_regions_for_homer.txt')
        
        subset = refseq[refseq['balb_nod_notx_1h_fc'] <= -1]
        yzer.run_homer(subset, 'promoter_overlap_notx_1h_nod_down', dirpath,
                       cpus=6, center=False, reverse=False, preceding=True, size=400, length=[8, 10, 12, 15], bg=bg)
        subset = refseq[refseq['balb_nod_notx_1h_fc'] >= 1]
        yzer.run_homer(subset, 'promoter_overlap_notx_1h_nod_up', dirpath,
                       cpus=6, center=False, reverse=False, preceding=True, size=400, length=[8, 10, 12, 15], bg=bg)
    # TSS
    if False:
        refseq = data[data['has_refseq'] == 1]
        refseq = refseq[refseq['transcript_score'] >= 4]
        if True:
            yzer.run_homer(refseq, 'refseq_tss', dirpath,
                       cpus=6, center=False, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

        bg = yzer.get_filename(dirpath, 'refseq_tss/refseq_tss_regions_for_homer.txt')
        
        subset = refseq[refseq['balb_nod_notx_1h_fc'] <= -1]
        yzer.run_homer(subset, 'tss_notx_1h_nod_down', dirpath,
                       cpus=6, center=False, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
        subset = refseq[refseq['balb_nod_notx_1h_fc'] >= 1]
        yzer.run_homer(subset, 'tss_notx_1h_nod_up', dirpath,
                       cpus=6, center=False, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
    # Enhancers
    if False:
        enh = data[data['distal'] == 't']
        enh = enh[enh['h3k4me2_tag_count'] > 10]
        enh = enh[enh[['nod_notx_1h_tag_count','balb_notx_1h_tag_count']].max(axis=1) >= 10]
        if False:
            yzer.run_homer(enh, 'enhancer_at_least_10_tags_15', dirpath,
                       cpus=6, center=False, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

        bg = yzer.get_filename(dirpath, 'enhancer_at_least_10_tags_15/enhancer_at_least_10_tags_15_regions_for_homer.txt')
        
        subset = enh[enh['balb_nod_notx_1h_fc'] <= -1]
        yzer.run_homer(subset, 'enhancer_at_least_10_tags_notx_1h_nod_down_15', dirpath,
                       cpus=6, center=False, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
        
        subset = enh[enh['balb_nod_notx_1h_fc'] >= 1]
        yzer.run_homer(subset, 'enhancer_at_least_10_tags_notx_1h_nod_up_15', dirpath,
                       cpus=6, center=False, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)

    # me2
    if True:
        enh = data[data['distal'] == 't']
        enh = enh[enh['h3k4me2_tag_count'] > 10]
        enh = enh[enh[['nod_notx_1h_tag_count','balb_notx_1h_tag_count']].max(axis=1) >= 10]
        peaks = yzer.import_file(yzer.get_filename(base_dirpath, 'transcripts_me2_peaks.txt'))
        enh = enh.merge(peaks, how='left', on='id')
        
        if True:
            yzer.run_homer(enh, 'me2_peaks', dirpath,
                       cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

        bg = yzer.get_filename(dirpath, 'me2_peaks/me2_peaks_regions_for_homer.txt')
        
        subset = enh[enh['balb_nod_notx_1h_fc'] <= -1]
        subset = subset.groupby('peak_id', as_index=False).mean()
        subset['transcription_start'] = subset['start']
        subset['transcription_end'] = subset['end']
        yzer.run_homer(subset, 'me2_peaks_notx_1h_nod_down_15', dirpath,
                       cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
        
        subset = enh[enh['balb_nod_notx_1h_fc'] >= 1]
        subset = subset.groupby('peak_id', as_index=False).mean()
        subset['transcription_start'] = subset['start']
        subset['transcription_end'] = subset['end']
        yzer.run_homer(subset, 'me2_peaks_notx_1h_nod_up_15', dirpath,
                       cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
