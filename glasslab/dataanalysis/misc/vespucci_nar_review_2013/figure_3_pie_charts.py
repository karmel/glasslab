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
    img_dirpath = yzer.get_and_create_path(dirpath, 'figure_3_pie_charts')
    
    yzer.legend_location = 'lower left'
    pie1 = '''Annotated by RefSeq and/or ncRNA.org    14,342
Unannotated    29,313'''
    pie1 = [row.split('    ') for row in pie1.split('\n')]
    pie1 = zip(*pie1)
    yzer.piechart(map(lambda s: int(s.replace(',','')), pie1[1]), 
                  pie1[0], 
                  title='Transcripts with Score >= 1', 
                  save_dir=img_dirpath, 
                  show_plot=True)
    
    pie2 = '''Promoter-associated RNA    5,623
Antisense of RefSeq    3,924
Post-TTS, same-strand    6,584
Other RefSeq Proximal    2,809
Distal with H3K4me1    6,393
Distal within 2kbp of H3K4me1    961
Remainder    3,019'''
    pie2 = [row.split('    ') for row in pie2.split('\n')]
    pie2 = zip(*pie2)
    yzer.legend_columns = 2
    yzer.piechart(map(lambda s: int(s.replace(',','')), pie2[1]), 
                  pie2[0], 
                  title='Unannotated Transcripts with Score >= 1', 
                  save_dir=img_dirpath, 
                  show_plot=True)