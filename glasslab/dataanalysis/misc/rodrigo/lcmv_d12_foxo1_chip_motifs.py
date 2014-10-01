'''
Created on Sep 27, 2014

@author: karmel
'''
from glasslab.dataanalysis.misc.rodrigo.samples import sample_name,\
    get_threshold
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer
if __name__ == '__main__':
    yzer = MotifAnalyzer()

    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/' +\
        'Miscellaneous_Collaborations/Rodrigo_CD8s_2014_09/Enhancers'
    dirpath = yzer.get_path(dirpath)

    cond, seq, breed = ('naive', 'foxo1', '')
    sample_prefix = sample_name(cond, seq, breed)
    sample_dirpath = yzer.get_filename(dirpath, sample_prefix)
    filename = yzer.get_filename(sample_dirpath,
                                 sample_prefix + '_enhancers.txt')

    data = yzer.import_file(filename)
    data = data.fillna(0)

    min_thresh = get_threshold(seq)
    data = data[data['tag_count'] >= min_thresh]

    fold = 2
    if True:
        naive_only = data[
            data['lcmv_d12_foxo1_tag_count'] < min_thresh]
        yzer.run_homer(naive_only,
                       'naive_only', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('naive_only', len(naive_only))

        shared = data[
            (data['tag_count'] * fold >= data['lcmv_d12_foxo1_tag_count']) &
            (data['lcmv_d12_foxo1_tag_count'] * fold >= data['tag_count'])]
        yzer.run_homer(shared,
                       'shared', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('shared', len(shared))

    cond, seq, breed = ('lcmv_d12', 'foxo1', '')
    sample_prefix = sample_name(cond, seq, breed)
    sample_dirpath = yzer.get_filename(dirpath, sample_prefix)
    filename = yzer.get_filename(sample_dirpath,
                                 sample_prefix + '_enhancers.txt')

    data = yzer.import_file(filename)
    data = data.fillna(0)

    min_thresh = get_threshold(seq)
    data = data[data['tag_count'] >= min_thresh]

    fold = 2
    if True:
        lcmv_d12_only = data[
            data['naive_foxo1_tag_count'] < min_thresh]
        yzer.run_homer(lcmv_d12_only,
                       'lcmv_d12_only', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('lcmv_d12_only', len(lcmv_d12_only))
