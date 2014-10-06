'''
Created on Oct 5, 2014

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
    breed_sets = get_breed_sets()

    # Working with only WT for now
    sample_names, short_names = breed_sets[0]

    # Set up our samples
    datasets = {}
    for sample_prefix in short_names:
        sample_dirpath = yzer.get_filename(dirpath, sample_prefix)
        filename = yzer.get_filename(sample_dirpath,
                                     sample_prefix + '_enhancers.txt')

        data = yzer.import_file(filename)
        data = data.fillna(0)

        min_thresh = get_threshold('atac')

        data = data[data['tag_count'] >= min_thresh]

        datasets[sample_prefix] = data

    # Now set up case-specific counts.

    # First, how many naive?
    print('Naive total:', len(datasets['d0']))
    print('Naive not d7:', sum(datasets['d0']['klrghi_d7']))
