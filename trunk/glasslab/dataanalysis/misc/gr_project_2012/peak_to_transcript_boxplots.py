'''
Created on Sep 7, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    grapher = SeqGrapher()
    base_dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    base_dirpath = grapher.get_path(base_dirpath)
    dirpath = grapher.get_filename(base_dirpath, 'motifs')
    filename = grapher.get_filename(dirpath, 'transcript_vectors.txt')
    
    data = grapher.import_file(filename)
    
    
    # Boxplots for gr_dex peaks by lfc in Dex
    if False:
        #data = data[data['distal'] == 't']
        data = data[data['has_refseq'] == 1]
        
        down = data[data['dex_1_lfc'] <= -1]
        up = data[data['dex_1_lfc'] >= 1]
        nc = data[abs(data['dex_1_lfc']) < 1]
        
        key = 'p65_kla_tag_count'
        datasets = [down[key],nc[key],up[key]]
        datasets = [d['p65_kla_dex_tag_count'] - d[key] for d in [down, nc, up]]
        
        #title = 'Tags in p65 peaks in KLA 1h + DMSO 2h: Distal'
        title = 'Diff in tags in p65 peaks in KLA 1h + Dex 2h vs KLA 1h + DMSO 2h: RefSeq'
        ax = grapher.boxplot(datasets, 
                        ['Down in Dex 2h','No change in Dex 2h','Up in Dex 2h',],
                         title=title, 
                         xlabel='Condition', 
                         ylabel='Total tags in all peaks overlapping transcript', 
                         show_outliers=False, show_plot=False)
        grapher.save_plot(grapher.get_filename(base_dirpath, 'boxplots', 'dex',
                               title.replace(' ','_')))
        grapher.show_plot()
        for sub in datasets: print sub.mean()
    
    # Boxplots for gr_kla_dex peaks by lfc in Dex
    if False:
        #data = data[data['distal'] == 't']
        #data = data[data['has_refseq'] == 1]
        
        trans = data[(data['kla_1_lfc'] >= 1) & (data['dex_over_kla_1_lfc'] <= -.58)]
        rest = data[(data['kla_1_lfc'] < 1) | (data['dex_over_kla_1_lfc'] > -.58)]
        
        key = 'gr_dex_tag_count'
        datasets = [rest[key],trans[key]]
        datasets = [d['gr_kla_dex_tag_count'] - d[key] for d in [rest, trans]]
        
        title = 'Tags in p65 peaks in KLA 1h + Dex 2h: Distal'
        title = 'Diff in tags in GR peaks in KLA 1h + Dex 2h vs Dex 1h'
        ax = grapher.boxplot(datasets, 
                        ['Not transrepressed in KLA 1h + Dex 2h','Transrepressed in KLA 1h + Dex 2h',],
                         title=title, 
                         xlabel='Condition', 
                         ylabel='Total tags in all peaks overlapping transcript', 
                         show_outliers=False, show_plot=False)
        grapher.save_plot(grapher.get_filename(base_dirpath, 'boxplots', 'kla_dex',
                               title.replace(' ','_')))
        grapher.show_plot()
        for sub in datasets: print sub.mean()
        
    if True:
        #data = data[data['has_refseq'] == 1]
        data = data[data['distal'] == 't']
        data['gr_diff'] = data['gr_kla_dex_tag_count'] - data['gr_dex_tag_count']
        data['p65_diff'] = data['p65_kla_dex_tag_count'] - data['p65_kla_tag_count']
        data['gr_by_length'] = data['gr_kla_dex_tag_count']/data['length']*10000
        data['p65_by_length'] = data['p65_kla_dex_tag_count']/data['length']*10000
        grapher.scatterplot(data, 'gr_kla_dex_tag_count', 'p65_diff',log=True)
        
        
        