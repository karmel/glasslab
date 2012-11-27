'''
Created on Oct 8, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from matplotlib import pyplot
from glasslab.dataanalysis.misc.gr_project_2012.boxplots_redistribution_pairs import get_high_quality_pairs

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

        data = data.groupby(['id','chr_name'],as_index=False).mean()
        
        gr_dex_gt_kla_dex = sum(data['tag_count_3']*1.5 < data['tag_count_4'])
        gr_dex_kla_gt_kla = sum(data['tag_count_4']*1.5 < data['tag_count_3'])
        has_gr = sum((data['tag_count_3'] != 0) | (data['tag_count_4'] != 0))
        no_gr = sum((data['tag_count_3'] == 0) & (data['tag_count_4'] == 0))
        
        counts = [gr_dex_gt_kla_dex, gr_dex_kla_gt_kla, (has_gr - gr_dex_gt_kla_dex - gr_dex_kla_gt_kla), no_gr]
        labels = ['Has more GR in Dex\nthan KLA+Dex', 'Has more GR in KLA+Dex\nthan Dex',
                  'No change in GR\nfrom Dex to KLA+Dex', 'No immediately overlapping GR']
        
        colors = ('#E5FCC2','#9DE0AD','#45ADA8','#547980')
        
        pyplot.figure(figsize=[10,10])
        patches, texts, autotexts = pyplot.pie(counts, labels=labels, colors=colors,
                                               autopct='%.2f%%', labeldistance=1.05)
        
        title = 'Distribution of GR at peaks that have 2x p65 in KLA+Dex vs. KLA'\
                + '\n(Total Regions = {0})'.format(len(data))
        pyplot.title(title)
        #pyplot.legend(loc='lower left')
        pyplot.savefig(yzer.get_filename(yzer.get_filename(dirpath, 'redistribution','distribution_of_gr_pie_chart_limited.png')))
        pyplot.show()
        