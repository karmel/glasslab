'''
Created on Oct 1, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_non_refseq_by_p65')
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'transcript_vectors.txt'))
    
    if True:
        data = transcripts[transcripts['refseq'] == 'f']
        data = data[data['has_infrastructure'] == 0]
        data = data[data['length'] < 6000]
        data = data[data['dex_1_lfc'] < 1]
        data = data[data['kla_1_lfc'] >= 1]
        data = data[data['gr_kla_dex_tag_count'] > 0]
        
        
        data = data.fillna(0)
        
        ratio = 3
        colname = 'dex_over_kla_1_lfc'
        none = (data['p65_kla_tag_count'] + data['p65_kla_dex_tag_count'] == 0)
        kla_gt = (data['p65_kla_tag_count'] > ratio*data['p65_kla_dex_tag_count'])
        kla_dex_gt = (data['p65_kla_dex_tag_count'] > ratio*data['p65_kla_tag_count'])
        nc = (data['p65_kla_tag_count'] + data['p65_kla_dex_tag_count'] > 0) \
            & (data['p65_kla_dex_tag_count'] < ratio*data['p65_kla_tag_count']) \
            & (data['p65_kla_tag_count'] < ratio*data['p65_kla_dex_tag_count'])
        
        
        data[kla_gt].to_csv(yzer.get_filename(img_dirpath, 'enhancer_like_lose_p65.txt'), 
                                      sep='\t', header=True, index=False)
        raise Exception
        title = 'LFC in KLA + Dex over KLA by change in p65:\nEnhancer-Like, Has GR in KLA+Dex'
        names = [s.format('p65') for s in ['No {0}','Loses {0} in KLA+Dex','No change in {0}', 'Gains {0} in KLA+Dex']]
        ax = yzer.boxplot([data[none][colname], data[kla_gt][colname], data[nc][colname], data[kla_dex_gt][colname]], 
                     names,
                     title=title, 
                     xlabel='p65 Status', 
                     ylabel='log2(KLA+Dex GRO-seq/KLA GRO-seq)', 
                     show_outliers=False, show_plot=False)
        yzer.save_plot(yzer.get_filename(img_dirpath, 'dex_over_kla_1_lfc_by_p65_sums_3x_change.png'))
        yzer.show_plot()