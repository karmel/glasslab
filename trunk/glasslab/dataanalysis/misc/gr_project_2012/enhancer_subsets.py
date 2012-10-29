'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
import numpy


if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'piecharts_by_binding')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_gene.txt'))
    
    # Make sure we have dimethyl
    #data = data[data.filter(like='h3k4me2').sum(axis=1) > 0]
    
    #data = yzer.collapse_strands(data)
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'feature_vectors.txt'))
    transcripts['id'] = transcripts['glass_transcript_id']
    data = data.merge(transcripts, how='left', on='id', suffixes=['','_trans'])
    
    data = data.fillna(0)
    
    transrepressed = data[(data['kla_1_lfc'] >= 1) & (data['dex_over_kla_1_lfc'] <= -.58)]
    not_trans = data[(data['kla_1_lfc'] < 1) | (data['dex_over_kla_1_lfc'] > -.58)]
    
    supersets = (('Transrepressed',transrepressed),('Not transrepressed', not_trans))
    
    # Plot trans versus not
    yzer.piechart([len(d) for d in zip(*supersets)[1]], zip(*supersets)[0],
                 title='Enhancer-like Subsets by state in KLA+Dex', 
                 save_dir=img_dirpath)
    
    tfs = [('PU.1','pu_1'),('p65','p65'),('GR','gr')]
    contexts = [('DMSO',''),('Dex','dex'),('KLA','kla'),('KLA+Dex','kla_dex')]
    
    for name, dataset in supersets:
        total_for_set = len(dataset)
        for tf_name, tf in tfs:
            # Get count for enhancer elements with this TF at all
            with_tf = dataset[dataset.filter(like=tf).sum(axis=1) > 0]
            without_tf = dataset[dataset.filter(like=tf).sum(axis=1) == 0]
            
            # Plot with TF versus not
            yzer.piechart([len(with_tf), len(without_tf)], 
                          ['Has {0}'.format(tf_name), 'No {0}'.format(tf_name)],
                         title='{0} Enhancer-like Subsets by {1} Occupancy in any state'.format(name, tf_name), 
                         save_dir=img_dirpath)
            
            
            
            