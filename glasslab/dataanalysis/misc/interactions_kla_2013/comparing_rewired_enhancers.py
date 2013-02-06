'''
Created on Jan 30, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from collections import OrderedDict

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_rewiring_lfc')
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'enhancer_sets','transcript_vectors.txt'))
    
    sets = OrderedDict((
             ('all', yzer.import_file(yzer.get_filename(data_dirpath,'all_vectors.cdt'))),
             #('all_6h', yzer.import_file(yzer.get_filename(data_dirpath,'kla_6h','all_vectors.cdt'))),
             ('rewired', yzer.import_file(yzer.get_filename(data_dirpath,'rewired_vectors.cdt'))),
             #('rewired_6h', yzer.import_file(yzer.get_filename(data_dirpath,'kla_6h','rewired_vectors.cdt'))),
             ('shared', yzer.import_file(yzer.get_filename(data_dirpath,'shared_vectors.cdt'))),
             ))

    for key, val in sets.items(): 
        sets[key + '_changed'] = sets[key][sets[key]['enhancer_lfc'].abs() > 2]
        sets[key + '_up'] = sets[key][sets[key]['enhancer_lfc'] > 2]
        sets[key + '_down'] = sets[key][sets[key]['enhancer_lfc'] < -2]
    
    if True:
        # Correlations
        
        for key, val in sets.iteritems(): 
            print key, len(val)
            print val.corr()
        
        notx_only = yzer.import_file(yzer.get_filename(data_dirpath,'notx_only_vectors.cdt'))
        kla_only = yzer.import_file(yzer.get_filename(data_dirpath,'kla_only_vectors.cdt'))
        print len(notx_only), len(kla_only)
        print sum(notx_only['enhancer_lfc'] < -1), sum(kla_only['enhancer_lfc'] < -1)
        print sum(notx_only['enhancer_lfc'] > 1), sum(kla_only['enhancer_lfc'] > 1)
        
        print sum(notx_only['notx_gene_lfc'] < -1), sum(kla_only['notx_gene_lfc'] < -1)
        print sum(notx_only['notx_gene_lfc'] > 1), sum(kla_only['notx_gene_lfc'] > 1)
        print sum(notx_only['kla_4h_gene_lfc'] < -1), sum(kla_only['kla_4h_gene_lfc'] < -1)
        print sum(notx_only['kla_4h_gene_lfc'] > 1), sum(kla_only['kla_4h_gene_lfc'] > 1)
    
    if True:
        # p65 occupancy
        for key, val in sets.iteritems(): 
            data = val.merge(transcripts, how='left', on='id')
            data = data.fillna(0)
            print key
            print data['h3k4me2_tag_count'].describe()
    
    


