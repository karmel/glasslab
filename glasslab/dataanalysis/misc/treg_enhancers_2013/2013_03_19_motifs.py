'''
Created on Feb 12, 2013

@author: karmel
'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    data_dirpath = 'karmel/Desktop/Projects/GlassLab/Data/Sequencing/ChipSeq/CD4TCell/2013_03_19'
    motif_dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/TReg_enhancers/motifs/2013_03_19'
    data_dirpath = yzer.get_path(data_dirpath)
    motif_dirpath = yzer.get_path(motif_dirpath)
    
    
    for cell_type in ('Naive','TReg'):
        for antibody in ('H3K4me2','H3K27Ac'):
            filename = yzer.get_filename(data_dirpath, 
                                         'FoxP3-GFP-CD4TCell-{0}-{1}-MNase-13-03-19'.format(cell_type,
                                                                                            antibody),
                                         'regions.txt')
            data = yzer.import_file(filename, skiprows=39)
            data = data.fillna(0)
            data['id'] = data['#PeakID']
            data['chr_name'] = data['chr']
            data['strand'] = (data['strand'] == '-').apply(int)

            if True:
                yzer.run_homer(data, 'all_{0}_{1}'.format(cell_type,antibody), motif_dirpath,
                           cpus=6, center=True, reverse=False, preceding=False, size=200, length=[8, 10, 12, 15])
    