'''
Created on Oct 8, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
import math

def get_high_quality_pairs(data, transcripts):
    data = data.merge(transcripts, how='left', on='glass_transcript_id',suffixes=['','trans'])
    data = data.fillna(0)
    
    
    
    data = data[(data['tag_count'] > 2*data['tag_count_2']) & # KLA+Dex p65 > KLA p65
                (data['reverse_tag_count'] > 2*data['reverse_tag_count_2'])]  # Nearby KLA p65 > KLA+Dex p65

    # Get rid of the cases where the two peaks are too different in size
    data = data[data.apply(lambda row: math.log(row['tag_count']/row['reverse_tag_count'],2), axis=1) <= 1.5]
    
    # Get rid of the cases where the transcript is up in Dex + KLA
    data = data[data['dex_1_lfc'] < 1]
    
    
    data['region_start'] = data.apply(lambda row: int(min(row['transcription_start'], row['transcription_start_5'])), axis=1)
    data['region_end'] = data.apply(lambda row: int(max(row['transcription_end'], row['transcription_end_5'])), axis=1)
    # Get rid of pairs that are really just overlapping
    data = data[data['region_end'] - data['region_start'] >= 300]
    #data = data[data['region_end'] - data['region_start'] <= 10000]
    
    return data 

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath,'motifs','from_peaks')
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'transcript_vectors.txt'))
    transcripts['glass_transcript_id'] = transcripts['id']
    
    if True:
        all_data = yzer.import_file(yzer.get_filename(dirpath, 'redistribution',
                                               'p65_peaks_bigger_in_kla_dex_with_nearby_bigger_kla_peaks.txt'))
        
        data = get_high_quality_pairs(all_data, transcripts)
        
        '''
        # Print these out to send to collaborators.
        data = data.groupby(['id','chr_name'],as_index=False).min()
        data['for_UCSC'] = data.apply(lambda row: '{0}:{1}-{2}'.format(row['chr_name'], 
                                                                       row['region_start'], 
                                                                       row['region_end']), axis=1)
        
        data.to_csv(yzer.get_filename(dirpath, 'redistribution', 'redistribution_peak_pairs_2x.txt'), 
                                      sep='\t', cols=['id', 'chr_name','region_start','region_end',
                                                      'for_UCSC', 'tag_count','tag_count_2',
                                                      'reverse_tag_count','reverse_tag_count_2'], 
                                      header=True, index=False)
        '''
        
        # Group data, take non-refseq transcripts preferably, and transrepressed preferably.
        data = data.sort(columns='dex_over_kla_1_lfc', axis=0)
        grouped = data.groupby(['id','chr_name'])

        def get_enhancer(group):
            enhancer = group[group['has_refseq'] == 1][:1]
            if not enhancer: enhancer = group[:1]
            return enhancer
        
        regrouped = grouped.apply(get_enhancer)
        regrouped = regrouped.reset_index(drop=True)
        
        
        regrouped = regrouped[regrouped['kla_1_lfc'] >=1]
        regrouped = regrouped[regrouped['has_refseq'] == 1]
        transcripts = transcripts[transcripts['has_refseq'] == 1]
        transcripts = transcripts[transcripts['kla_1_lfc'] >= 1]
        
        all_trans = transcripts['dex_over_kla_1_lfc']
        with_p65 = transcripts[transcripts['p65_kla_tag_count'] > 0]['dex_over_kla_1_lfc']
        with_pair = regrouped['dex_over_kla_1_lfc']
        title = 'Transcript log fold change at redistribution pairs:\nRefseq, up in KLA alone'
        names = ['All transcripts', 'Transcripts with p65', 
                 'Transcripts with\na redistribution pair']
        ax = yzer.boxplot([all_trans, with_p65, with_pair], 
                     names,
                     title=title, 
                     xlabel='Transcript subset', 
                     ylabel='log2(KLA+Dex / KLA)', 
                     show_outliers=False, show_plot=False)
        yzer.save_plot(yzer.get_filename(dirpath, 'redistribution', 'boxplots', 
                                         'dex_over_kla_lfc_boxplot_1_limited_refseq_up_in_kla.png'))
        yzer.show_plot()