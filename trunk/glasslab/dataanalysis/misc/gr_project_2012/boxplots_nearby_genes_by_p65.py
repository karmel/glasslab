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
    img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_nearby_genes_by_p65')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'transcript_vectors.txt'))
    nearby = yzer.import_file(yzer.get_filename(img_dirpath, 'nearest_genes_to_enhancer_like_less_p65.txt'))
    
    if True:
        data = data.fillna(0)
        nearby = nearby.fillna(0)
        
        data = data[data['has_refseq'] == 1]
        data = data[~data['id'].isin(nearby['id'])]
        
        nearby = nearby.groupby(['id'],as_index=False).mean()
        ratio = 3
        colname = 'dex_over_kla_1_lfc'
        none = (data['p65_kla_tag_count'] + data['p65_kla_dex_tag_count'] == 0)
        kla_gt = (data['p65_kla_tag_count'] > ratio*data['p65_kla_dex_tag_count'])
        kla_dex_gt = (data['p65_kla_dex_tag_count'] > ratio*data['p65_kla_tag_count'])
        nc = (data['p65_kla_tag_count'] + data['p65_kla_dex_tag_count'] > 0) \
            & (data['p65_kla_dex_tag_count'] < ratio*data['p65_kla_tag_count']) \
            & (data['p65_kla_tag_count'] < ratio*data['p65_kla_dex_tag_count'])
        
        
        #data[kla_gt].to_csv(yzer.get_filename(img_dirpath, 'enhancer_like_lose_p65.txt'), 
        #                              sep='\t', header=True, index=False)
        
        title = 'LFC in KLA + Dex over KLA by change in p65:\nRefSeq'
        names = [s.format('p65') for s in ['No {0}','Loses {0}\nin KLA+Dex','No change in {0}', 'Gains {0}\nin KLA+Dex',
                                           'Near Enhancer\nthat loses {0}']]
        
        groups = [data[none][colname], data[kla_gt][colname], 
                           data[nc][colname], data[kla_dex_gt][colname],
                           nearby[colname]]
        for g in groups: print len(g)
        ax = yzer.boxplot(groups, 
                     names,
                     title=title, 
                     xlabel='Transcript Status', 
                     ylabel='log2(KLA+Dex GRO-seq/KLA GRO-seq)', 
                     show_outliers=False, show_plot=False)
        yzer.save_plot(yzer.get_filename(img_dirpath, 'dex_over_kla_1_lfc_with_nearby_unique_3x_change.png'))
        yzer.show_plot()