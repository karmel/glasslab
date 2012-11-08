'''
Created on Nov 7, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.dataanalysis.misc.cd4tcell_finland_2012.resources import replicate_sets

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells_Finland_2012'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'basic_scatterplots')

    data = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
                            
    for rep_set in replicate_sets:
        (key1a, key2a), (key1b, key2b) = rep_set
        
        # Make sure we get at least 10 tags in one of each replicate pair
        dataset = data[((data['{0}_tag_count'.format(key1a)] > 10) | (data['{0}_tag_count'.format(key2a)] > 10))
                       & ((data['{0}_tag_count'.format(key1b)] > 10) | (data['{0}_tag_count'.format(key2b)] > 10))]
        
        
        for lfc in (1,2):
            min_lfc = lfc
            
            lfc_key_a = '{0}_{1}_lfc'.format(key1a, key2a)
            lfc_key_b = '{0}_{1}_lfc'.format(key1b, key2b)
            # Make sure one set is above thresh, the other is at least at min_lfc
            up = dataset[((dataset[lfc_key_a] >= lfc) & (dataset[lfc_key_b] >= min_lfc))
                         | ((dataset[lfc_key_b] >= lfc) & (dataset[lfc_key_a] >= min_lfc))]
            down = dataset[((dataset[lfc_key_a] <= -lfc) & (dataset[lfc_key_b] <= -min_lfc))
                         | ((dataset[lfc_key_b] <= -lfc) & (dataset[lfc_key_a] <= -min_lfc))]
            
            print '\n\n\nUp {0}x in {1} and {2} versus {3} and {4}:'.format(2**lfc, 
                                                                    key2a, key2b, key1a, key1b)
            print ','.join(up['gene_names'])
            print '\n\n\nDown {0}x in {1} and {2} versus {3} and {4}:'.format(2**lfc, 
                                                                    key2a, key2b, key1a, key1b)
            print ','.join(down['gene_names'])
        