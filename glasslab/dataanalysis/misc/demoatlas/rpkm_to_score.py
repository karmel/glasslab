'''
Created on Feb 5, 2013

@author: karmel

Note: Made font.weight = bold and axes.titlesize = 24, font.size = 16 in matplotlibrc
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from matplotlib.ticker import FormatStrFormatter
from matplotlib import pyplot

class PrettyAxisGrapher(SeqGrapher):
    def other_plot(self, ax=None, *args, **kwargs):
        major_formatter = FormatStrFormatter('%d')
        ax.xaxis.set_major_formatter(major_formatter)
        ax.yaxis.set_major_formatter(major_formatter)
        ax.xaxis.set_ticks([2**x for x in (1, 5, 10, 15, 20)])
        ax.yaxis.set_ticks([2**x for x in (1, 5, 10, 15, 20)])
        return ax
    
    def show_count_scatterplot(self, data, ax, text_shift=0, text_color='black'):
        pyplot.text(.125, .825+text_shift, 'Total count: %d' % len(data),
                        color=text_color,
                        transform=ax.transAxes)
        
    def show_correlation_scatterplot(self, data, xcolname, ycolname, ax, text_shift=0, text_color='black'):
        pyplot.text(.125, .8+text_shift, 'r = %.4f' % data[xcolname].corr(data[ycolname]),
                    color=text_color,
                    transform=ax.transAxes)
        
if __name__ == '__main__':
    yzer = PrettyAxisGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/Demo-data'
    dirpath = yzer.get_path(dirpath)
    
    img_dirpath = yzer.get_and_create_path(dirpath, 'scatterplots/large_gap_500bp')
    
    data = yzer.import_file(yzer.get_filename(dirpath,'score_rpkm_500bp.txt'))
    ax = yzer.scatterplot(data, 
                 xcolname='rpkm', ycolname='score', log=True,  
                 title='Transcript RPKM versus Vespucci Score', 
                 xlabel='RPKM',
                 ylabel='Vespucci Transcript Score',  
                 show_2x_range=False, plot_regression=False, 
                 show_count=True, show_correlation=True, show_legend=False, 
                 save_dir=img_dirpath, show_plot=True)
    ax = yzer.scatterplot(data, 
                 xcolname='length', ycolname='score', log=True,  
                 title='Transcript Length versus Vespucci Score', 
                 xlabel='Length (bp)',
                 ylabel='Vespucci Transcript Score',  
                 show_2x_range=False, plot_regression=False, 
                 show_count=True, show_correlation=True, show_legend=False, 
                 save_dir=img_dirpath, show_plot=True)
    ax = yzer.scatterplot(data, 
                 xcolname='length', ycolname='rpkm', log=True,  
                 title='Transcript Length versus RPKM', 
                 xlabel='Length (bp)',
                 ylabel='RPKM',  
                 show_2x_range=False, plot_regression=False, 
                 show_count=True, show_correlation=True, show_legend=False, 
                 save_dir=img_dirpath, show_plot=True)