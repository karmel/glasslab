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
             ('all_6h', yzer.import_file(yzer.get_filename(data_dirpath,'kla_6h','all_vectors.cdt'))),
             ('rewired', yzer.import_file(yzer.get_filename(data_dirpath,'rewired_vectors.cdt')))
             ('rewired_6h', yzer.import_file(yzer.get_filename(data_dirpath,'kla_6h','rewired_vectors.cdt')))))

    sets['rewired_changed'] = sets['rewired'][sets['rewired']['enhancer_lfc'].abs() > 1]
    sets['rewired_up'] = sets['rewired'][sets['rewired']['enhancer_lfc'] > 1]
    sets['rewired_down'] = sets['rewired'][sets['rewired']['enhancer_lfc'] < -1]
    
    if True:
        # Correlations
        
        for key, val in sets.iteritems(): 
            print key
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
            print key
            print data['p65_tag_count'].mean(), data['p65_tag_count'].median(), data['p65_tag_count'].quantiles(4)
    
    


