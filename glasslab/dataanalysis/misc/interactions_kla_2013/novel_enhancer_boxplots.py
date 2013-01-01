'''
Created on Dec 31, 2012

@author: karmel

We want to compare the number of interactions going in to enhancers 
that are new in KLA (as determined by H3K4me2 being absent in notx)
versus those that are already existant. 

Here then we separate transcripts that are distal and have me2
according to whether they have less than a fraction of the average me2 in notx
versus all of the KLA runs (1h, 6h, 24h, 48h), or whether they have
at least the average me2.

We count the interactions connected to each transcript and
draw a boxplot. In order to most easily pull in all the 
zero-interaction enhancers, we load those with a separate query
and use the difference in count.
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/novel_enhancers'
    dirpath = yzer.get_path(dirpath)
    
    img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_by_me2_in_notx')

    data = {}
    data['all_enhancers_less_me2'] = yzer.import_file(yzer.get_filename(dirpath, 'all_enhancers_less_me2_in_notx.txt'))
    data['all_enhancers_with_me2'] = yzer.import_file(yzer.get_filename(dirpath, 'all_enhancers_with_me2_in_notx.txt'))
    data['interacting_in_notx_less_me2'] = yzer.import_file(yzer.get_filename(dirpath, 'interacting_in_notx_enhancers_less_me2_in_notx.txt'))
    data['interacting_in_notx_with_me2'] = yzer.import_file(yzer.get_filename(dirpath, 'interacting_in_notx_enhancers_with_me2_in_notx.txt'))
    data['interacting_in_kla_less_me2'] = yzer.import_file(yzer.get_filename(dirpath, 'interacting_in_kla_enhancers_less_me2_in_notx.txt'))
    data['interacting_in_kla_with_me2'] = yzer.import_file(yzer.get_filename(dirpath, 'interacting_in_kla_enhancers_with_me2_in_notx.txt'))
    data['interacting_in_kla_4h_less_me2'] = yzer.import_file(yzer.get_filename(dirpath, 'interacting_in_kla_4h_enhancers_less_me2_in_notx.txt'))
    data['interacting_in_kla_4h_with_me2'] = yzer.import_file(yzer.get_filename(dirpath, 'interacting_in_kla_4h_enhancers_with_me2_in_notx.txt'))
   
    for key in data: 
        if key[:5] != 'inter' or 'less_me2' not in key: continue
        data[key] = data[key][data[key]['me2_notx_tag_count']*4 < data[key]['me2_tag_count']]
    
    zero_intxns_in_notx_less_me2 = len(data['all_enhancers_less_me2']) - len(data['interacting_in_notx_less_me2'])
    zero_intxns_in_notx_with_me2 = len(data['all_enhancers_with_me2']) - len(data['interacting_in_notx_with_me2'])
    zero_intxns_in_kla_less_me2 = len(data['all_enhancers_less_me2']) - len(data['interacting_in_kla_less_me2'])
    zero_intxns_in_kla_with_me2 = len(data['all_enhancers_with_me2']) - len(data['interacting_in_kla_with_me2'])
    zero_intxns_in_kla_4h_less_me2 = len(data['all_enhancers_less_me2']) - len(data['interacting_in_kla_4h_less_me2'])
    zero_intxns_in_kla_4h_with_me2 = len(data['all_enhancers_with_me2']) - len(data['interacting_in_kla_4h_with_me2'])
    
    notx_less_me2_counts = data['interacting_in_notx_less_me2']['count'].values.tolist() + [0]*zero_intxns_in_notx_less_me2
    notx_with_me2_counts = data['interacting_in_notx_with_me2']['count'].values.tolist() + [0]*zero_intxns_in_notx_with_me2
    kla_less_me2_counts = data['interacting_in_kla_less_me2']['count'].values.tolist() + [0]*zero_intxns_in_kla_less_me2
    kla_with_me2_counts = data['interacting_in_kla_with_me2']['count'].values.tolist() + [0]*zero_intxns_in_kla_with_me2
    kla_4h_less_me2_counts = data['interacting_in_kla_4h_less_me2']['count'].values.tolist() + [0]*zero_intxns_in_kla_4h_less_me2
    kla_4h_with_me2_counts = data['interacting_in_kla_4h_with_me2']['count'].values.tolist() + [0]*zero_intxns_in_kla_4h_with_me2

    labels = ['Less than 1/4\navg H3K4me2 in notx,\ninteractions in notx', 
              'At least avg H3K4me2 in notx\ninteractions in notx',
              'Less than 1/4\navg H3K4me2 in notx,\ninteractions in KLA 30m', 
              'At least avg H3K4me2 in notx\ninteractions in KLA 30m',
              'Less than 1/4\navg H3K4me2 in notx,\ninteractions in KLA 4h', 
              'At least avg H3K4me2 in notx\ninteractions in KLA 4h',]
    vals = [notx_less_me2_counts, notx_with_me2_counts, 
            kla_less_me2_counts, kla_with_me2_counts,
            kla_4h_less_me2_counts, kla_4h_with_me2_counts]
    labels = [l + '\n(count: {0})'.format(len(v)) for l,v in zip(labels, vals)]
    
    title = 'Number of interactions with "enhancers" by H3K4me2 state in notx'
    ax = yzer.boxplot(vals, labels, 
                     title=title, xlabel='Enhancer subset', 
                     ylabel='Number of interactions with other transcripts', 
                     show_outliers=False, show_plot=True, wide=True,
                     save_dir=img_dirpath)
