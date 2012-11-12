'''
Created on Nov 7, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.misc.cd4tcell_finland_2012.resources import replicate_sets
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer

if __name__ == '__main__':
    yzer = MotifAnalyzer()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells_Finland_2012'
    dirpath = yzer.get_path(dirpath)
    curr_path = yzer.get_and_create_path(dirpath, 'motif_analysis')

    data = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    
    if False:
        yzer.run_homer(data, 'all_refseq_preceding', curr_path, 
                center=False, reverse=False, preceding=True, size=200,
                cpus=6)
        yzer.run_homer(data, 'all_refseq', curr_path, 
                center=False, reverse=False, preceding=False, size=200,
                cpus=6)
    bg_preceding = yzer.get_filename(curr_path, 'all_refseq_preceding','all_refseq_preceding_regions_for_homer.txt')
    bg = yzer.get_filename(curr_path, 'all_refseq','all_refseq_regions_for_homer.txt')
               
    for rep_set in replicate_sets:
        (key1a, key2a), (key1b, key2b) = rep_set
        
        # Make sure we get at least 10 tags in one of each replicate pair
        dataset = data[((data['{0}_tag_count'.format(key1a)] > 10) | (data['{0}_tag_count'.format(key2a)] > 10))
                       & ((data['{0}_tag_count'.format(key1b)] > 10) | (data['{0}_tag_count'.format(key2b)] > 10))]
        
        
        
        for lfc in (1,):
            min_lfc = lfc/2
            
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
        
            yzer.run_homer(up, 'up_{0}x_{1}_{2}_preceding'.format(2**lfc, lfc_key_a, lfc_key_b), curr_path, 
                    center=False, reverse=False, preceding=True, size=200,
                    cpus=6, bg=bg_preceding)
            yzer.run_homer(up, 'up_{0}x_{1}_{2}'.format(2**lfc, lfc_key_a, lfc_key_b), curr_path, 
                    center=False, reverse=False, preceding=False, size=200,
                    cpus=6, bg=bg)
            yzer.run_homer(up, 'up_{0}x_{1}_{2}'.format(2**lfc, lfc_key_a, lfc_key_b), curr_path, 
                    center=True, reverse=False, preceding=False, size=300,
                    cpus=6, bg=bg, files_already_prepped=True)
            yzer.run_homer(down, 'down_{0}x_{1}_{2}_preceding'.format(2**lfc, lfc_key_a, lfc_key_b), curr_path, 
                    center=False, reverse=False, preceding=True, size=200,
                    cpus=6, bg=bg_preceding)
            yzer.run_homer(down, 'down_{0}x_{1}_{2}'.format(2**lfc, lfc_key_a, lfc_key_b), curr_path, 
                    center=False, reverse=False, preceding=False, size=200,
                    cpus=6, bg=bg)
            yzer.run_homer(down, 'down_{0}x_{1}_{2}'.format(2**lfc, lfc_key_a, lfc_key_b), curr_path, 
                    center=True, reverse=False, preceding=False, size=300,
                    cpus=6, bg=bg, files_already_prepped=True)
