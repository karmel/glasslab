'''
Created on Oct 1, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
import random

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_nearby_genes_by_p65')
    
    
    ratio = 1.5
        
    data = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'transcript_vectors.txt'))
    nearby = yzer.import_file(yzer.get_filename(img_dirpath, 'nearest_genes_to_enhancer_like_less_p65_{0}x.txt'.format(ratio)))
    colname = 'dex_over_kla_1_lfc'
        
    
    pausing = True
    if pausing:
        colname = 'pausing_ratio_ratio'
        
        # We want previously calculated bucket scores,
        # Joined to old transcripts because we have since updated IDs
        bucket_scores = yzer.import_file(yzer.get_filename(dirpath, 'classification','feature_vectors.txt'))
        old_transcripts = yzer.import_file(yzer.get_filename(dirpath, 'classification','transcript_vectors_old_ids.txt'))
        bucket_scores['id'] = bucket_scores['glass_transcript_id']
        bucket_scores[colname] = bucket_scores['kla_dex_1_bucket_score']/bucket_scores['kla_1_bucket_score']
        
        bucket_scores = bucket_scores.merge(old_transcripts, how='left', on='id',suffixes=['','_old'])
        bucket_scores = bucket_scores[['id','gene_names','pausing_ratio_ratio']]
        
        # Clean up gene names
        key = 'gene_names_clean'
        for dataset in (bucket_scores, nearby, data):
            dataset[key] = dataset['gene_names'].map(lambda x: str(x).replace('{','').replace('}','')) 
            
        nearby = nearby.merge(bucket_scores, how='left', on=key,suffixes=['','_bucket'])
        data = data.merge(bucket_scores, how='left', on=key,suffixes=['','_bucket'])
        
    if True:
        data = data.fillna(0)
        nearby = nearby.fillna(0)
        
        data = data[data['has_refseq'] == 1]
        data = data[~data['id'].isin(nearby['id'])]
        
        nearby = nearby.groupby(['id'],as_index=False).mean()
        
        none = (data['p65_kla_tag_count'] + data['p65_kla_dex_tag_count'] == 0)
        kla_gt = (data['p65_kla_tag_count'] > ratio*data['p65_kla_dex_tag_count'])
        kla_dex_gt = (data['p65_kla_dex_tag_count'] > ratio*data['p65_kla_tag_count'])
        nc = (data['p65_kla_tag_count'] + data['p65_kla_dex_tag_count'] > 0) \
            & (data['p65_kla_dex_tag_count'] < ratio*data['p65_kla_tag_count']) \
            & (data['p65_kla_tag_count'] < ratio*data['p65_kla_dex_tag_count'])
        
        
        names = [s.format('p65') for s in ['No {0}','Loses {0}\nin KLA+Dex','No change in {0}', 'Gains {0}\nin KLA+Dex',
                                           'Near Enhancer\nthat loses {{0}} ({0}-fold)'.format(ratio)]]
        
        groups = [data[none], data[kla_gt], data[nc], data[kla_dex_gt]]
        
        # We want to randomly sample to get equi-sized groups
        desired = len(nearby)
        for i, g in enumerate(groups):
            rows = random.sample(g.index, desired) 
            groups[i] = g.ix[rows]
            
        to_plot = [g[colname] for g in (groups + [nearby])]
        
        
        title = 'LFC in KLA + Dex over KLA by change in p65:' \
                    + '\nRefSeq, randomly sampled to {0} transcripts'.format(desired)
        
        if pausing: 
            title = 'Pausing Ratio Ratio by change in p65:' \
                        + '\nRefSeq, randomly sampled to {0} transcripts'.format(desired)
        
        ax = yzer.boxplot(to_plot, 
                     names,
                     title=title, 
                     xlabel='Transcript Status', 
                     ylabel=(pausing and 'PausingRatio(KLA+Dex)/PausingRatio(KLA)')\
                        or 'log2(KLA+Dex GRO-seq/KLA GRO-seq)', 
                     show_outliers=False, show_plot=False)
        yzer.save_plot(yzer.get_filename(img_dirpath, 
                '{2}_with_nearby_unique_{0}x_change_sampled_{1}.png'.format(ratio, random.randint(0,9999),
                                                                            colname)))
        yzer.show_plot()