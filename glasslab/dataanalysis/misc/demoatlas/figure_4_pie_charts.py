'''
Created on Apr 19, 2013

@author: karmel

Note: Made font.weight = bold and line.width = 4 in matplotlibrc
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/Demo-data'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'figure_4_pie_charts')
    
    pie1 = '''Annotated by RefSeq and/or ncRNA.org    17,882
Unannotated    46,354'''
    pie1 = [row.split('    ') for row in pie1.split('\n')]
    pie1 = zip(*pie1)
    yzer.piechart(map(lambda s: int(s.replace(',','')), pie1[1]), 
                  pie1[0], 
                  title='Transcripts with Score >= 2', 
                  save_dir=img_dirpath, 
                  show_plot=True)
    
    pie2 = '''Refseq proximal    23,600
Distal with H3K4me1    7,670
Distal within 2kbp of H3K4me1    3,802
Within 5kbp of RefSeq TTS    5,206
Remainder    6,076'''
    pie2 = [row.split('    ') for row in pie2.split('\n')]
    pie2 = zip(*pie2)
    yzer.piechart(map(lambda s: int(s.replace(',','')), pie2[1]), 
                  pie2[0], 
                  title='Unannotated Transcripts with Score >= 2', 
                  save_dir=img_dirpath, 
                  show_plot=True)