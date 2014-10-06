'''
Created on Oct 6, 2014

@author: karmel
'''
from glasslab.dataanalysis.misc.rodrigo.samples import get_threshold,\
    get_breed_sets
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer
if __name__ == '__main__':
    yzer = MotifAnalyzer()

    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/' +\
        'Miscellaneous_Collaborations/Rodrigo_CD8s_2014_09/Enhancers_set2'
    dirpath = yzer.get_path(dirpath)

    datasets = {}
    samples, short_names = get_breed_sets()[0]

    for j, sample_prefix in enumerate(short_names):
        sample_dirpath = yzer.get_filename(dirpath, sample_prefix)
        filename = yzer.get_filename(sample_dirpath,
                                     sample_prefix + '_enhancers.txt')

        data = yzer.import_file(filename)
        data = data.fillna(0)

        min_thresh = get_threshold('atac')

        data = data[data['tag_count'] >= min_thresh]

        datasets[sample_prefix] = data

    naive_ct = len(datasets['d0'])
    for celltype in ('klrghi', 'klrglo'):
        d7_with_naive = sum(
            datasets[celltype + '_d7']['d0_tag_count'] >= min_thresh)
        d7_no_naive = sum(
            datasets[celltype + '_d7']['d0_tag_count'] < min_thresh)

        d12_with_d7 = sum(
            datasets[celltype + '_d12'][celltype + '_d7_tag_count'] >= min_thresh)
        d12_no_d7 = datasets[celltype + '_d12'][
            datasets[celltype + '_d12'][celltype + '_d7_tag_count'] < min_thresh]
        d12_naive = sum(d12_no_d7['d0_tag_count'] >= min_thresh)
        d12_no_naive = sum(d12_no_d7['d0_tag_count'] < min_thresh)

        d35_with_d7 = sum(
            datasets[celltype + '_d35'][celltype + '_d7_tag_count'] >= min_thresh)
        d35_no_d7 = datasets[celltype + '_d35'][
            datasets[celltype + '_d35'][celltype + '_d7_tag_count'] < min_thresh]
        d35_naive = sum(d35_no_d7['d0_tag_count'] >= min_thresh)
        d35_no_naive = sum(d35_no_d7['d0_tag_count'] < min_thresh)

        print(naive_ct)
        print(d7_with_naive, d7_no_naive)
        print(d12_with_d7, d12_naive, d12_no_naive)
        print(d35_with_d7, d35_naive, d35_no_naive)
