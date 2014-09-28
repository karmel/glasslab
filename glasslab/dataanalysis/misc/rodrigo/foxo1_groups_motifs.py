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

    cond, seq, breed = ('naive', 'atac', '')
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
        # ATAC peaks that are absent in the FOXO1 KO
        foxo1_critical = data[
            data['foxo1_ko_naive_atac_tag_count'] < min_thresh]
        yzer.run_homer(foxo1_critical,
                       'foxo1_critical', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('foxo1_critical', len(foxo1_critical))

        # ATAC peaks exist, but are less than half as big
        foxo1_important = data[
            (data['foxo1_ko_naive_atac_tag_count'] >= min_thresh) &
            (data['foxo1_ko_naive_atac_tag_count'] * fold < data['tag_count'])]
        yzer.run_homer(foxo1_important,
                       'foxo1_important', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('foxo1_important', len(foxo1_important))

        # The sum of foxo1_critical and foxo1_important
        foxo1_dependent = data[
            (data['foxo1_ko_naive_atac_tag_count'] * fold < data['tag_count'])]
        yzer.run_homer(foxo1_dependent,
                       'foxo1_dependent', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('foxo1_dependent', len(foxo1_dependent))

        # ATAC peaks that don't change with KO of Foxo1
        foxo1_independent = data[
            (data['foxo1_ko_naive_atac_tag_count'] * fold >= data['tag_count'])]
        yzer.run_homer(foxo1_independent,
                       'foxo1_independent', sample_dirpath,
                       cpus=10, center=True, reverse=False, preceding=False,
                       size=200, length=[8, 10, 12], mock=True)
        print('foxo1_independent', len(foxo1_independent))
