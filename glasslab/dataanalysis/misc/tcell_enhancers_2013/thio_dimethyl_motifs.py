'''
Created on Feb 12, 2013

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/H3K4me2/Analysis'
    dirpath = yzer.get_path(dirpath)
    motif_dirpath = yzer.get_filename(dirpath, 'motifs')
    filename = yzer.get_filename(dirpath, 'thio_peak_vectors.txt')
    data = yzer.import_file(filename)
    data = data.fillna(0)


    # me2
    if True:
        if False:
            yzer.run_homer(data, 'thio_all', motif_dirpath,
                       cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])
        
        data = data[data['tss_id'] == 0]
        if True:
            yzer.run_homer(data, 'thio_h3k4me2_distal', motif_dirpath,
                       cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])
        
        subset = data[data['ctcf_tag_count'] > 0]
        if False:
            yzer.run_homer(subset, 'thio_h3k4me2_distal_with_ctcf', motif_dirpath,
                       cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])
        
 