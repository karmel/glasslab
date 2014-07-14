'''
Created on Jul 2, 2014

@author: karmel

We want to compare acetylation at various subsets of me2 enhancers.
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/AND_TCR/Analysis/Chips1_2_3/Motifs'
    dirpath = yzer.get_path(dirpath)
    
    savepath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/AND_TCR/Analysis/Chips1_2_3/Figures/Acetyl_box_plots'
    savepath = yzer.get_path(savepath)

    if False:
        peptide, ab = 'PCC', 'me2'
        pep_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(peptide, ab))
            
        filename = yzer.get_filename(pep_dirpath, 
                        '{}_{}_enhancers_with_ac.txt'.format(peptide, ab))

        data = yzer.import_file(filename)
        data = data.fillna(0)
        
        de_novo = data[data['no_pep_tag_count'] == 0]
        established = data[data['no_pep_tag_count'] > 0]
        established_batf = established[established['batf_tag_count'] > 0]
        de_novo_batf = de_novo[de_novo['batf_tag_count'] > 0]
        yzer.boxplot([established['pcc_ac_tag_count'], 
                      de_novo['pcc_ac_tag_count'],
                      established_batf['pcc_ac_tag_count'],
                      de_novo_batf['pcc_ac_tag_count'],],
                     bar_names=['Established', 'De novo',
                                'Established with Batf',
                                'De novo with Batf'])
        
        
    if False:
        peptide, ab = 'K99A', 'me2'
        pep_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(peptide, ab))
            
        filename = yzer.get_filename(pep_dirpath, 
                        '{}_{}_enhancers_with_ac.txt'.format(peptide, ab))

        data = yzer.import_file(filename)
        data = data.fillna(0)
        
        de_novo = data[data['no_pep_tag_count'] == 0]
        established = data[data['no_pep_tag_count'] > 0]
        established_rel = established[established['rel_id'] > 0]
        de_novo_rel = de_novo[de_novo['rel_id'] > 0]
        yzer.boxplot([established['k99a_ac_tag_count'], 
                      de_novo['k99a_ac_tag_count'],
                      established_rel['k99a_ac_tag_count'],
                      de_novo_rel['k99a_ac_tag_count'],],
                     bar_names=['Established', 'De novo',
                                'Established with REL motif',
                                'De novo with REL motif'])
    if False:
        peptide, ab = 'NoPep', 'me2'
        pep_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(peptide, ab))
            
        filename = yzer.get_filename(pep_dirpath, 
                        '{}_{}_enhancers_with_ac.txt'.format(peptide, ab))

        data = yzer.import_file(filename)
        data = data.fillna(0)
        
        de_novo = data[(data['k99a_tag_count'] == 0) & (data['pcc_tag_count'] == 0)]
        established = data[(data['k99a_tag_count'] > 0) | (data['pcc_tag_count'] > 0)]
        established_rel = established[established['rel_id'] > 0]
        de_novo_rel = de_novo[de_novo['rel_id'] > 0]
        yzer.boxplot([established['no_pep_ac_tag_count'], 
                      de_novo['no_pep_ac_tag_count'],
                      established_rel['no_pep_ac_tag_count'],
                      de_novo_rel['no_pep_ac_tag_count'],],
                     bar_names=['Established', 'De novo',
                                'Established with REL motif',
                                'De novo with REL motif'])
        
    if True:
        peptide, ab = 'NoPep', 'ac'
        pep_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(peptide, ab))
            
        filename = yzer.get_filename(pep_dirpath, 
                        '{}_{}_enhancers_with_rel.txt'.format(peptide, ab))

        data = yzer.import_file(filename)
        data = data.fillna(0)
        
        with_rel = data[data['rel_id'] > 0]
        
        yzer.boxplot([with_rel['tag_count'], 
                      with_rel['k99a_tag_count'],
                      with_rel['pcc_tag_count'],
                      ],
                     bar_names=['NoPep H3K27Ac', 'K99A H3K27Ac',
                                'PCC H3K27Ac',],
                     xlabel='', ylabel='Tag count',
                     title='H3K27Ac Tag Counts at REL motifs in CD4-Null EnNoPephancers')