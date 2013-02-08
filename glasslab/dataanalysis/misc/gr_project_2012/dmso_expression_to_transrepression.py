'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from matplotlib import pyplot
from glasslab.utils.functions import nonzero
from glasslab.dataanalysis.misc.gr_project_2012.v1.elongation import total_tags_per_run

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    
    consistent = False
    img_dirpath = yzer.get_and_create_path(dirpath, 'dmso_expression_to_transrepression', consistent and 'consistent' or 'rep1')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_gene.txt'))
    data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
    
    draw_pies = True
    min_tags = 30
    ratio = 1.5
    # Make sure we have dimethyl
    data = data[data.filter(like='h3k4me2').max(axis=1) > min_tags]
    data = data[data['minimal_distance'] >= 1000]
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    transcripts['nearest_refseq_transcript_id'] = transcripts['id']
    data = data.merge(transcripts, how='left', on='nearest_refseq_transcript_id', suffixes=['','_trans'])
    
    data = data.fillna(0)
    
    total_tags = total_tags_per_run()
    data['dmso_1_rpkm'] = data['dmso_1_tag_count_trans']*(10**3*10**6)/data['length_trans']/total_tags['dmso'][1]
    data['h4k8ac_kla_ratio'] = nonzero(data['h4k8ac_kla_tag_count'])/nonzero(data['h4k8ac_tag_count'])
    data['h4k8ac_kla_dex_ratio'] = nonzero(data['h4k8ac_kla_dex_tag_count'])/nonzero(data['h4k8ac_kla_tag_count'])
    
    for subgroup, suffix, dataset in (('RefSeq Transcripts', '_trans', 
                                       data.groupby(by='nearest_refseq_transcript_id', as_index=False).mean()),):
        
        ax = yzer.scatterplot(dataset[(dataset['kla_1_lfc_trans'] >= 1)], 
                              'dmso_1_rpkm', 
                              'dex_over_kla_1_lfc_trans', log=True, 
                              title='GR transrepression by DMSO expression for Up-regulated genes', 
                              xlabel='DMSO 2h RPKM', 
                              ylabel='log2(KLA+Dex GRO-seq/DMSO GRO-seq)',
                              show_2x_range=False, plot_regression=True, show_count=True, 
                              show_correlation=True, save_dir=img_dirpath, show_plot=True)
        ax = yzer.scatterplot(dataset[(dataset['kla_1_lfc_trans'] >= 1)], 
                              'h4k8ac_kla_ratio', 
                              'dex_over_kla_1_lfc_trans', log=True, 
                              title='GR transrepression by KLA to DMSO H4K8ac tag ratio for Up-regulated genes', 
                              xlabel='KLA Tags/DMSO Tags', 
                              ylabel='log2(KLA+Dex GRO-seq/DMSO GRO-seq)',
                              show_2x_range=False, plot_regression=True, show_count=True, 
                              show_correlation=True, save_dir=img_dirpath, show_plot=True)
        ax = yzer.scatterplot(dataset[(dataset['kla_1_lfc_trans'] >= 1)], 
                              'h4k8ac_kla_dex_ratio', 
                              'dex_over_kla_1_lfc_trans', log=True, 
                              title='GR transrepression by KLA+Dex to KLA H4K8ac tag ratio for Up-regulated genes', 
                              xlabel='Dex+KLA Tags/KLA Tags', 
                              ylabel='log2(KLA+Dex GRO-seq/DMSO GRO-seq)',
                              show_2x_range=False, plot_regression=True, show_count=True, 
                              show_correlation=True, save_dir=img_dirpath, show_plot=True)