'''
Created on Feb 12, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/NOD_BALBc/ThioMacs/Analysis_2013_02/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'srf_binding')
    data = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    
    data = data.fillna(0)
    data = data[data[['nod_notx_1h_tag_count','balb_notx_1h_tag_count']].max(axis=1) >= 10]

    subsets = [data, data[(data['has_refseq'] == 1) & (data['transcript_score'] >= 4)], 
                   data[(data['distal'] == 't') & (data['h3k4me2_tag_count'] > 10)]]    
    
    
    # Add in nearest genes for enhancers
    enh = subsets[2].copy()
    nearest_genes = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_genes.txt'))
    enh['enh_id'] = enh['id']
    enh = enh.merge(nearest_genes, how='left', on='enh_id', suffixes=['_enh',''])
    # Now id refers to the nearest gene; use that to create the vector of balb_nod_notx_1h_fc
    # Keep the original srf_notx_tag_count though, to refer to enhancer
    gene_ids = enh[['id','srf_notx_tag_count']]
    gene_ids = gene_ids.merge(data, how='left', on='id', suffixes=['','_gene'])
    
    subsets.append(gene_ids)
    
    datasets = []
    for subset in subsets:
        with_srf = subset[subset['srf_notx_tag_count'] > 10]
        without_srf = subset[subset['srf_notx_tag_count'] == 0]
        datasets += [with_srf, without_srf]
        
    
    
    vals = [d['balb_nod_notx_1h_fc'] for d in datasets]
    labels = ['All with SRF', 'All without SRF', 'RefSeq with SRF', 'RefSeq without SRF', 
              'Enhancers with SRF', 'Enhancers without SRF', 
              'Genes Nearest to\nEnhancers with SRF','Genes Nearest to\nEnhancers without SRF']
    labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(labels, vals)]
    title = 'Balbc vs NOD Log-fold Change According to SRF Binding'
    if True:
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel='Transcript subset', 
                         ylabel='log2(NOD notx 1h GRO-seq/BALBc notx 1h GRO-seq)', 
                         show_outliers=False, show_plot=True, wide=True,
                         save_dir=img_dirpath)
    if False:
        yzer = MotifAnalyzer()
        motifs_dirpath = yzer.get_filename(dirpath, 'motifs')
        # Motif finding for enhancers with SRF
        bg = yzer.get_filename(motifs_dirpath, 'enhancer_at_least_10_tags_15/enhancer_at_least_10_tags_15_regions_for_homer.txt')
        
        enh_with_srf = datasets[4]
        yzer.run_homer(enh_with_srf, 'enhancer_with_srf', motifs_dirpath,
                       cpus=6, center=False, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
    
        # Motif finding for enhancers with SRF and down
        yzer.run_homer(enh_with_srf[enh_with_srf['balb_nod_notx_1h_fc'] <= -1], 'enhancer_with_srf_down', motifs_dirpath,
                       cpus=6, center=False, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
    
        # Motif finding for enhancers with SRF and down versus SRF
        bg = yzer.get_filename(motifs_dirpath, 'enhancer_with_srf/enhancer_with_srf_regions_for_homer.txt')
        yzer.run_homer(enh_with_srf[enh_with_srf['balb_nod_notx_1h_fc'] <= -1], 'enhancer_with_srf_down_vs_with_srf', motifs_dirpath,
                       cpus=6, center=False, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
