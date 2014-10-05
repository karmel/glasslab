'''
Created on Sep 27, 2014

@author: karmel
'''
from glasslab.dataanalysis.misc.rodrigo.samples import sample_name,\
    get_threshold, get_breed_sets
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer
if __name__ == '__main__':
    yzer = MotifAnalyzer()

    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/' +\
        'Miscellaneous_Collaborations/Rodrigo_CD8s_2014_09/Enhancers_set2'
    dirpath = yzer.get_path(dirpath)

    for samples, short_names in get_breed_sets():
        for sample_prefix in short_names:
            sample_dirpath = yzer.get_filename(dirpath, sample_prefix)
            filename = yzer.get_filename(sample_dirpath,
                                         sample_prefix + '_enhancers.txt')

            data = yzer.import_file(filename)
            data = data.fillna(0)

            if True:

                min_thresh = get_threshold('atac')

                subdata = data[data['tag_count'] >= min_thresh]
                yzer.run_homer(subdata,
                               'all', sample_dirpath,
                               cpus=10, center=True, reverse=False, preceding=False,
                               size=200, length=[8, 10, 12], mock=True)
