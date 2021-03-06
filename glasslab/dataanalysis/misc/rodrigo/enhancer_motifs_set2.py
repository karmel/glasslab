'''
Created on Sep 27, 2014

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

            if False:

                yzer.run_homer(data,
                               'all', sample_dirpath,
                               cpus=10, center=True, reverse=False, preceding=False,
                               size=200, length=[8, 10, 12], mock=True)
            if True:
                # Versus comparable sample in other breed
                subdata = data[
                    data['{}_tag_count'.format(oth_breed[1][j])] < min_thresh]
                yzer.run_homer(subdata,
                               'not_in_' + oth_breed[1][j], sample_dirpath,
                               cpus=10, center=True, reverse=False,
                               preceding=False,
                               size=200, length=[8, 10, 12], mock=True)

                # Versus hi/lo sample and prior sample if not naive.
                if j > 0:
                    if 'klrghi' in sample_prefix:
                        other = sample_prefix.replace('hi', 'lo')
                    elif 'klrglo' in sample_prefix:
                        other = sample_prefix.replace('lo', 'hi')
                    subdata = data[
                        data['{}_tag_count'.format(other)] < min_thresh]
                    yzer.run_homer(subdata,
                                   'not_in_' +
                                   other, sample_dirpath,
                                   cpus=10, center=True, reverse=False,
                                   preceding=False,
                                   size=200, length=[8, 10, 12], mock=True)

                    # One further back in timecourse
                    other = ('d7' in sample_prefix and short_names[0])\
                        or short_names[j - 2]  # Previous timepoint, skip hi/lo
                    subdata = data[
                        data['{}_tag_count'.format(other)] < min_thresh]
                    yzer.run_homer(subdata,
                                   'not_in_' +
                                   other, sample_dirpath,
                                   cpus=10, center=True, reverse=False,
                                   preceding=False,
                                   size=200, length=[8, 10, 12], mock=True)
