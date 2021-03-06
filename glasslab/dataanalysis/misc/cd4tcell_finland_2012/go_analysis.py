'''
Created on Nov 7, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.misc.cd4tcell_finland_2012.resources import replicate_sets
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer

if __name__ == '__main__':
    yzer = MotifAnalyzer()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells_Finland_2012/Analysis_2013_02'
    dirpath = yzer.get_path(dirpath)
    go_path = yzer.get_and_create_path(dirpath, 'with_me3','go_analysis', '0_8_min_lfc')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    data = data.fillna(0)
    data = data[data['naive_me3_tag_count'] + data['act_me3_tag_count'] > 0]
    
    if False:
        curr_path = yzer.get_and_create_path(dirpath, 'with_me3', 'motif_analysis')

        yzer.run_homer(data, 'all_refseq_preceding', curr_path, 
                center=False, reverse=False, preceding=True, size=200,
                cpus=6)
        yzer.run_homer(data, 'all_refseq', curr_path, 
                center=False, reverse=False, preceding=False, size=200,
                cpus=6)
        yzer.run_homer(data, 'all_refseq', curr_path, 
                center=True, reverse=False, preceding=False, size=300,
                cpus=6, files_already_prepped=True)
        
        bg_preceding = yzer.get_filename(curr_path, 'all_refseq_preceding','all_refseq_preceding_regions_for_homer.txt')
        bg = yzer.get_filename(curr_path, 'all_refseq','all_refseq_regions_for_homer.txt')
    
    gene_file = file(yzer.get_filename(go_path, 'differential_gene_names.txt'), 'w')          
    for rep_set in replicate_sets:
        (key1a, key2a), (key1b, key2b) = rep_set
        
        # Make sure we get at least 10 tags in one of each replicate pair
        dataset = data[((data['{0}_tag_count'.format(key1a)] > 10) | (data['{0}_tag_count'.format(key2a)] > 10))
                       & ((data['{0}_tag_count'.format(key1b)] > 10) | (data['{0}_tag_count'.format(key2b)] > 10))]
        
        
        
        for lfc in (1,):
            min_lfc = .8
            
            lfc_key_a = '{0}_{1}_lfc'.format(key1a, key2a)
            lfc_key_b = '{0}_{1}_lfc'.format(key1b, key2b)
            # Make sure one set is above thresh, the other is at least at min_lfc
            up = dataset[((dataset[lfc_key_a] >= lfc) & (dataset[lfc_key_b] >= min_lfc))
                         | ((dataset[lfc_key_b] >= lfc) & (dataset[lfc_key_a] >= min_lfc))]
            down = dataset[((dataset[lfc_key_a] <= -lfc) & (dataset[lfc_key_b] <= -min_lfc))
                         | ((dataset[lfc_key_b] <= -lfc) & (dataset[lfc_key_a] <= -min_lfc))]
            
            gene_file.write('\n\n\nUp {0}x in {1} and {2} versus {3} and {4}:'.format(2**lfc, 
                                                                    key2a, key2b, key1a, key1b))
            gene_file.write(','.join(up['gene_names']))
            gene_file.write('\n\n\nDown {0}x in {1} and {2} versus {3} and {4}:'.format(2**lfc, 
                                                                    key2a, key2b, key1a, key1b))
            gene_file.write(','.join(down['gene_names']))
        
            if False:
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

    gene_file.close()