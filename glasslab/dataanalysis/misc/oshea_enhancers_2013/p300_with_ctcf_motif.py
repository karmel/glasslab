'''
Created on Feb 12, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/Oshea_enhancers/peak_overlaps'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'figures')
    
    peak = 'p300'
    th1 = yzer.import_file(yzer.get_filename(dirpath, 'th1_with_th2_{0}.txt'.format(peak))).fillna(0)
    th2 = yzer.import_file(yzer.get_filename(dirpath, 'th2_with_th1_{0}.txt'.format(peak))).fillna(0)
    th1_with_ctcf_motif = yzer.import_file(yzer.get_filename(dirpath, 'motifs','th_p300_enhancers_ctcf',
                                                         'th1_only_{0}_with_ctcf_motif.txt'.format(peak)))
    th2_with_ctcf_motif = yzer.import_file(yzer.get_filename(dirpath, 'motifs','th_p300_enhancers_ctcf',
                                                         'th2_only_{0}_with_ctcf_motif.txt'.format(peak)))
    shared_with_ctcf_motif = yzer.import_file(yzer.get_filename(dirpath, 'motifs','th_p300_enhancers_ctcf',
                                                         'th_shared_{0}_with_ctcf_motif.txt'.format(peak)))
    
    # Filter out promoters
    th1 = th1[th1['tss_id'] == 0]
    th2 = th2[th2['tss_id'] == 0]
    th1_with_ctcf_motif['id'] = th1_with_ctcf_motif['PositionID']
    th2_with_ctcf_motif['id'] = th2_with_ctcf_motif['PositionID']
    shared_with_ctcf_motif['id'] = shared_with_ctcf_motif['PositionID']
    
    th1['th1_tag_count'] = nonzero(th1['tag_count'])
    th1['th2_tag_count'] = nonzero(th1['p2_tag_count'])
    th2['th1_tag_count'] = nonzero(th2['tag_count'])
    th2['th2_tag_count'] = nonzero(th2['p2_tag_count'])
    

    with_ctcf = th1[th1['id'].isin(th1_with_ctcf_motif['id']) | th1['id'].isin(shared_with_ctcf_motif['id'])]
    without_ctcf = th1[~th1['id'].isin(th1_with_ctcf_motif['id']) & ~th1['id'].isin(shared_with_ctcf_motif['id'])]
    
    
    datasets = [with_ctcf, without_ctcf, th1, th1[th1['p2_tag_count'] > 0]]#, th2]
        
    
    vals = [d['th2_tag_count']/d['th1_tag_count'] for d in datasets]
    
    labels = ['With p300\nand CTCF motif', 'With p300\nbut not CTCF motif',
              'All in Th1', 'All Shared', 'All in Th2']
    
    labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(labels, vals)]
    title = 'Th1 and Th2 p300 Tag Count Ratios According to CTCF Motif Presence'
    if True:
        ax = yzer.boxplot(vals[:2], labels[:2], 
                         title=title, xlabel='Enhancer subset', 
                         ylabel='Tags in p300 peak in Th2/Tags in p300 peak in Th1', 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        ax = yzer.boxplot(vals, labels, 
                         title=title + '- All groups', xlabel='Enhancer subset', 
                         ylabel='Tags in p300 peak in Th2/Tags in p300 peak in Th1', 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)

        