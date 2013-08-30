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
    pie1 = '''Annotated by RefSeq and/or ncRNA.org    14,309
Unannotated    71,513'''
    pie1 = [row.split('    ') for row in pie1.split('\n')]
    pie1 = zip(*pie1)
    yzer.piechart(map(lambda s: int(s.replace(',','')), pie1[1]), 
                  pie1[0], 
                  title='MCF-7 Transcripts with Score >= 1', 
                  save_dir=img_dirpath, 
                  show_plot=True)
    
    pie2 = '''Promoter-associated RNA    7,145
Antisense of RefSeq    7,946
Post-TTS, same-strand    6,052
Other RefSeq Proximal    11,582
Distal with H3K4me2    2,352
Distal within 2kbp of H3K4me2    5,524
Distal remainder with LINE    16,292
Remainder    14,620'''
    pie2 = [row.split('    ') for row in pie2.split('\n')]
    pie2 = zip(*pie2)
    yzer.legend_columns = 2
    yzer.piechart(map(lambda s: int(s.replace(',','')), pie2[1]), 
                  pie2[0], 
                  title='Unannotated MCF-7 Transcripts with Score >= 1', 
                  save_dir=img_dirpath, 
                  show_plot=True)