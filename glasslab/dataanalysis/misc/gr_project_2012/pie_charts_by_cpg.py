'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'cpg_island_promoters', 'piecharts','by_genes_with_gr')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancer_classification','enhancers_with_nearest_gene.txt'))
    data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
    
    min_tags = 30
    # Make sure we have dimethyl
    data = data[data.filter(like='h3k4me2').max(axis=1) > min_tags]
    data = data[data['minimal_distance'] >= 1000]
    
    #data = yzer.collapse_strands(data)
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'cpg_island_promoters', 'transcript_vectors.txt'))
    transcripts['nearest_refseq_transcript_id'] = transcripts['id']
    data = data.merge(transcripts, how='left', on='nearest_refseq_transcript_id', suffixes=['','_trans'])
    data = data.merge(transcripts, how='left', on='id', suffixes=['','_enh'])
    
    data = data.fillna(0)
    
    transrepressed = data[(data['kla_1_lfc_trans'] >= 1) & (data['dex_over_kla_1_lfc_trans'] <= -.58)]
    not_trans = data[(data['kla_1_lfc_trans'] < 1) | (data['dex_over_kla_1_lfc_trans'] > -.58)]
    up_in_kla = data[(data['kla_1_lfc_trans'] >= 1) & (data['dex_over_kla_1_lfc_trans'] > -.58)]
    
    supersets = (('All', data),
                 ('Not near transrepressed genes', not_trans), 
                 ('Up in KLA', up_in_kla),
                 ('Near transrepressed genes',transrepressed))
    
    # Plot trans versus not
    yzer.piechart([len(d) for d in zip(*supersets[1:])[1]], zip(*supersets[1:])[0],
                 title='Enhancer-like Subsets by state in KLA+Dex', 
                 save_dir=img_dirpath, show_plot=False)
    
    for name, dataset in supersets:
        total_for_set = len(dataset)
        dataset = dataset[dataset['gr_kla_dex_tag_count'] > min_tags]
        # Get count for enhancer elements with/out CpG
        with_cpg = dataset[dataset['has_cpg_enh'] == 1]
        without_cpg = dataset[dataset['has_cpg_enh'] == 0]
        
        # Plot with TF versus not
        yzer.piechart([ len(without_cpg), len(with_cpg)], 
                      ['No CpG Island', 'Has CpG Island'],
                     title='Enhancer-like Subsets {0}\nby Overlap with CpG Island'.format(name.title()), 
                     save_dir=img_dirpath, show_plot=False)
            