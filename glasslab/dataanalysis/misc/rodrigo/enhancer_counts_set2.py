'''
Created on Sep 27, 2014

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.dataanalysis.misc.rodrigo.samples import get_threshold,\
    get_breed_sets
if __name__ == '__main__':
    yzer = SeqGrapher()

    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/' +\
        'Miscellaneous_Collaborations/Rodrigo_CD8s_2014_09/Enhancers_set2'
    dirpath = yzer.get_path(dirpath)

    save_path = yzer.get_and_create_path(
        dirpath, 'Figures', 'Enhancer_counts')

    datasets = {}
    breed_sets = get_breed_sets()
    for i, (samples, short_names) in enumerate(breed_sets):
        oth_breed = breed_sets[1 - i]
        for j, sample_prefix in enumerate(short_names):
            sample_dirpath = yzer.get_filename(dirpath, sample_prefix)
            filename = yzer.get_filename(sample_dirpath,
                                         sample_prefix + '_enhancers.txt')

            data = yzer.import_file(filename)
            data = data.fillna(0)

            min_thresh = get_threshold('atac')

            data = data[data['tag_count'] >= min_thresh]

            datasets[sample_prefix] = data

    # How many denovo d7 enhancers are also in foxo1 kos?
    for celltype in ('hi', 'lo'):
        d7 = datasets['klrg{}_d7'.format(celltype)]
        de_novo = d7[d7['d0_tag_count'] < min_thresh]

        all_shared = d7[
            'foxo1_ko_klrg{}_d7_tag_count'.format(celltype)] >= min_thresh
        all_not_shared = d7[
            'foxo1_ko_klrg{}_d7_tag_count'.format(celltype)] < min_thresh
        shared = de_novo[
            'foxo1_ko_klrg{}_d7_tag_count'.format(celltype)] >= min_thresh
        not_shared = de_novo[
            'foxo1_ko_klrg{}_d7_tag_count'.format(celltype)] < min_thresh

        labels = ['Also in Foxo1 KO', 'Not in Foxo1 KO']
        yzer.piechart([sum(all_shared), sum(all_not_shared)],
                      labels,
                      title='WT KLRG{} d7 Enhancers'.format(celltype),
                      save_dir=save_path, show_plot=False)
        yzer.piechart([sum(shared), sum(not_shared)],
                      labels,
                      title='WT KLRG{} d7 De Novo Enhancers'.format(celltype),
                      save_dir=save_path, show_plot=False)

        yzer.boxplot([d7[all_shared]['tag_count'].tolist(),
                      d7[all_not_shared]['tag_count'].tolist()],
                     labels,
                     title='ATAC-seq tags in WT KLRG{} d7 Enhancers'.format(
                         celltype),
                     ylabel='ATAC peak tag count', save_dir=save_path,
                     show_plot=False)
        yzer.boxplot([de_novo[shared]['tag_count'].tolist(),
                      de_novo[not_shared]['tag_count'].tolist()],
                     labels,
                     title='ATAC-seq tags in WT KLRG{} d7 De Novo Enhancers'.format(
                         celltype),
                     ylabel='ATAC peak tag count', save_dir=save_path,
                     show_plot=False)
