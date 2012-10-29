'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'piecharts_by_binding')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_gene.txt'))
    
    #data = yzer.collapse_strands(data)
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'feature_vectors.txt'))
    
    data = data.merge(transcripts, how='left outer', suffixes=['','_trans'])
    
    transrepressed = data[(data['kla_1_lfc'] >= 1) & (data['dex_over_kla_1_lfc'] <= -.58)]
    not_trans = data[(data['kla_1_lfc'] < 1) | (data['dex_over_kla_1_lfc'] > -.58)]
    
    total = len(data)
    
    tfs = [('PU.1','pu_1'),('p65','p65'),('GR','gr')]
    contexts = [('DMSO',''),('Dex','dex'),('KLA','kla'),('KLA+Dex','kla_dex')]
    
    for name, dataset in (('Transrepressed',transrepressed), 
                          ('Not transrepressed'), not_trans):
        total_for_set = len(dataset)
        for tf_name, tf in tfs:
            with_tf = 
            total_with_tf = len()
            
            
            