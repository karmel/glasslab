'''
Created on Feb 12, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/Oshea_enhancers/peak_overlaps'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'figures')
    
    peak_pretty = 'p300'
    peak = peak_pretty.lower()
    th1 = yzer.import_file(yzer.get_filename(dirpath, 'th1_with_th2_{0}.txt'.format(peak))).fillna(0)
    th2 = yzer.import_file(yzer.get_filename(dirpath, 'th2_with_th1_{0}.txt'.format(peak))).fillna(0)
    
    # Filter out promoters
    th1 = th1[th1['tss_id'] == 0]
    th2 = th2[th2['tss_id'] == 0]
    
    # Get venn-diagram sets
    only_th1 = th1[th1['p2_id'] == 0]
    only_th2 = th2[th2['p2_id'] == 0]
    shared = th1[th1['p2_id'] != 0]
    shared_check = th2[th2['p2_id'] != 0]
    print len(only_th1), len(only_th2), len(shared), len(shared_check)
    
    only_th1['th1_tag_count'] = nonzero(only_th1['tag_count'])
    only_th1['th2_tag_count'] = nonzero(only_th1['p2_tag_count'])
    shared['th1_tag_count'] = nonzero(shared['tag_count'])
    shared['th2_tag_count'] = nonzero(shared['p2_tag_count'])
    only_th2['th1_tag_count'] = nonzero(only_th2['p2_tag_count'])
    only_th2['th2_tag_count'] = nonzero(only_th2['tag_count'])
    
    data = shared.append(only_th1, ignore_index=True)
    data = data.append(only_th2, ignore_index=True)
    
    if False:
        # Scatterplots of tag counts
        ax = yzer.scatterplot(data, 'th1_tag_count', 'th2_tag_count',
                                log=True, color='blue', 
                                title='Th1 versus Th2 {0} Tag Counts at Peaks'.format(peak_pretty),
                                show_2x_range=True, show_legend=False, plot_regression=False,
                                show_count=True, show_correlation=True, 
                                save_dir=img_dirpath, show_plot=True)
    
    if True:
        # Motif finding.
        yzer = MotifAnalyzer()
        motifs_dirpath = yzer.get_and_create_path(dirpath, 'motifs')
        
        data['id'] = data.index
        
        if False:
            yzer.run_homer(data, 'th_{0}_enhancers'.format(peak), motifs_dirpath,
                       cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])

        bg = yzer.get_filename(motifs_dirpath, 'th_{0}_enhancers/th_{0}_enhancers_regions_for_homer.txt'.format(peak))
        multiple = 2
        
        '''
        yzer.run_homer(only_th1, 'th1_only_{0}_enhancers'.format(peak), motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
        yzer.run_homer(only_th2, 'th2_only_{0}_enhancers'.format(peak), motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
        yzer.run_homer(shared, 'th_shared_{0}_enhancers'.format(peak), motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
        
        yzer.run_homer(data[data['th1_tag_count'] > data['th2_tag_count']*multiple], 
                       'th1_up_{1}x_{0}_enhancers'.format(peak, multiple), motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
        yzer.run_homer(data[data['th2_tag_count'] > data['th1_tag_count']*multiple], 
                       'th2_up_{1}x_{0}_enhancers'.format(peak, multiple), motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
        '''
        yzer.run_homer(data[(data['th1_tag_count'] > data['th2_tag_count']*multiple) & (data['ctcf_tag_count'] > 0)], 
                       'th1_up_{1}x_{0}_enhancers_with_ctcf'.format(peak, multiple), motifs_dirpath,
                   cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15], bg=bg)
        
        
        
        