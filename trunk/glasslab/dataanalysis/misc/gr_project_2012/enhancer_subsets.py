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
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'feature_vectors.txt'))
    
    data = data.merge(transcripts, how='left outer', suffixes=['','_trans'])