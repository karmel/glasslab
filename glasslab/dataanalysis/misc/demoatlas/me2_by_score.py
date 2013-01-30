'''
Created on Jan 15, 2013

@author: karmel

Scatterplots of H3K4me2 peak tag counts by GROseq score
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/Demo-data'
    dirpath = yzer.get_path(dirpath)
    
    img_dirpath = yzer.get_and_create_path(dirpath, 'scatterplots')
    
    data = yzer.import_file(yzer.get_filename(dirpath,'me2_peaks_with_transcripts.txt'))
    data = data.fillna(0)
    data = data.groupby(by='id', as_index=True).mean()
    data['transcript_score'] = data['score(2)']
    ax = yzer.scatterplot(data, 
                 xcolname='transcript_score', ycolname='tag_count', log=True,  
                 title='H3K4me2 Tag Count as a Function of Transcript Score', 
                 xlabel='Glass Atlas Transcript Score', 
                 ylabel='Normalized H3K4me2 tag count', 
                 show_2x_range=True, plot_regression=True, 
                 show_count=True, show_correlation=True, show_legend=False, 
                 save_dir=img_dirpath, show_plot=True)
    