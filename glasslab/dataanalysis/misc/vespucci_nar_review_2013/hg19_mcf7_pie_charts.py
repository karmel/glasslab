'''
Created on Apr 19, 2013

@author: karmel

Note: Made font.weight = bold and axes.titlesize = 24 in matplotlibrc
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/NAR_review_data/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'hg19_mcf7_pie_charts')
    
    yzer.legend_location = 'lower left'
    pie1 = '''Annotated by RefSeq and/or ncRNA.org    14,022
Unannotated    67,046'''
    pie1 = [row.split('    ') for row in pie1.split('\n')]
    pie1 = zip(*pie1)
    yzer.piechart(map(lambda s: int(s.replace(',','')), pie1[1]), 
                  pie1[0], 
                  title='Hah et al MCF-7 Transcripts\nwith Score >= 1', 
                  save_dir=img_dirpath, 
                  show_plot=True)
    
    pie2 = '''Promoter-associated RNA    7,055
Antisense of RefSeq    7,539
Other RefSeq Proximal    13,664
Distal with H3K4me2    2,352
Distal w/in 2kbp of H3K4me2    5,524
Distal remainder with LINE    16,292
Remainder    14,620'''
    pie2 = [row.split('    ') for row in pie2.split('\n')]
    pie2 = zip(*pie2)
    yzer.legend_columns = 2
    yzer.piechart(map(lambda s: int(s.replace(',','')), pie2[1]), 
                  pie2[0], 
                  title='Hah et al Unannotated MCF-7 Transcripts\nwith Score >= 1', 
                  save_dir=img_dirpath, 
                  show_plot=True)