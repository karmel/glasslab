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
    
    # Filter out promoters
    th1 = th1[th1['tss_id'] == 0]
    th2 = th2[th2['tss_id'] == 0]
    
    
    th1['th1_tag_count'] = nonzero(th1['tag_count'])
    th1['th2_tag_count'] = nonzero(th1['p2_tag_count'])
    th2['th1_tag_count'] = nonzero(th2['tag_count'])
    th2['th2_tag_count'] = nonzero(th2['p2_tag_count'])
    

    with_ctcf = th1[th1['ctcf_tag_count'] > 0]
    without_ctcf = th1[th1['ctcf_tag_count'] == 0]
    
    
    datasets = [with_ctcf, without_ctcf, th1, th1[th1['p2_tag_count'] > 0]]
        
    
    vals = [d['th2_tag_count']/d['th1_tag_count'] for d in datasets]
    
    base_labels = ['With p300 and CTCF', 'With p300 but not CTCF',
                   'All in Th1', 'All Shared', 'All in Th2']
    
    labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(base_labels, vals)]
    title = 'Th1 and Th2 {0} Tag Count Ratios According to CTCF Co-localization'.format(peak)
    if False:
        ax = yzer.boxplot(vals[:2], labels[:2], 
                         title=title, xlabel='Enhancer subset', 
                         ylabel='Tags in {0} peak in Th2/Tags in {0} peak in Th1'.format(peak), 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
        ax = yzer.boxplot(vals, labels, 
                         title=title + '- All groups', xlabel='Enhancer subset', 
                         ylabel='Tags in {0} peak in Th2/Tags in {0} peak in Th1'.format(peak), 
                         show_outliers=False, show_plot=True, wide=False,
                         save_dir=img_dirpath)
    
    if True:
        datasets = [with_ctcf, without_ctcf, 
                   without_ctcf[without_ctcf['tag_count'] > 2*without_ctcf['p2_tag_count']],
                   th1, th1[th1['p2_tag_count'] > 0]]
        base_labels = ['With p300 and CTCF', 'With p300 but not CTCF',
                      'With p300 but not CTCF,\ndown in Th2',
                      'All in Th1', 'All Shared']

        for tf in ('Stat1KO', 'Stat4KO', 'TbetKO'):
            vals = [nonzero(d['{0}_{1}_tag_count'.format(tf.lower(), peak)])\
                    /nonzero(d['th1_tag_count']) for d in datasets]
            
            labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(base_labels, vals)]
            title = 'Th1 and {0} {1} Tag Count Ratios According to CTCF Co-localization'.format(tf, peak)
            ax = yzer.boxplot(vals[:2], labels[:2], 
                             title=title, xlabel='Enhancer subset', 
                             ylabel='Tags in {0} peak in {1}/Tags in {0} peak in Th1'.format(peak, tf), 
                             show_outliers=False, show_plot=False, wide=False,
                             save_dir=img_dirpath)
            ax = yzer.boxplot(vals, labels, 
                             title=title + '- All groups', xlabel='Enhancer subset', 
                             ylabel='Tags in {0} peak in {1}/Tags in {0} peak in Th1'.format(peak, tf), 
                             show_outliers=False, show_plot=False, wide=True,
                             save_dir=img_dirpath)
            
