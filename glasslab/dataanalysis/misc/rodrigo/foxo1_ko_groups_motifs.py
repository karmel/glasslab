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

    cond, seq, breed = ('naive', 'atac', 'foxo1_ko_')
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
        # ATAC peaks that are absent in the FOXO1 competent
        ko_new = data[
            data['naive_atac_tag_count'] < min_thresh]
        yzer.run_homer(ko_new,
                       'ko_new', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('ko_new', len(ko_new))

        # ATAC peaks exist, but are less than half as big
        ko_grows = data[
            (data['naive_atac_tag_count'] >= min_thresh) &
            (data['naive_atac_tag_count'] * fold < data['tag_count'])]
        yzer.run_homer(ko_grows,
                       'ko_grows', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('ko_grows', len(ko_grows))

        # The sum of new and bigger
        ko_dependent = data[
            (data['naive_atac_tag_count'] * fold < data['tag_count'])]
        yzer.run_homer(ko_dependent,
                       'ko_dependent', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('ko_dependent', len(ko_dependent))

        # ATAC peaks that don't change without KO of Foxo1
        ko_independent = data[
            (data['naive_atac_tag_count'] * fold >= data['tag_count'])]
        yzer.run_homer(ko_independent,
                       'ko_independent', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('ko_independent', len(ko_independent))
