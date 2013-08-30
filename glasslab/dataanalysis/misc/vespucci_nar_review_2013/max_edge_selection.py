'''
Created on Aug 27, 2013

@author: karmel

We want to evaluate various values of MAX_EDGE using the GRO5 data from 
Michael Lam's paper. We want as many transcripts as possible to contain
exactly one start site, and we want to penalize the merging of multiple start
sites. We calculate a success rate accordingly, using data from a number of 
different builds of Vespucci, and plot.

Note: Made font.weight = bold and axes.titlesize = 24  and axes.labelweight=bold in matplotlibrc
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from collections import OrderedDict
import pandas


if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/NAR_review_data/parameters'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'success_rates')
    
    edges = [100, 200, 500, 800, 1000, 1500, 2000, 5000]
    datasets = OrderedDict()
    for e in edges:
        datasets[e] = pandas.read_csv(yzer.get_filename(dirpath, 'vs_gro5',
                                                    'gro5_count_me_{}.txt'.format(e)),
                                      header=0, 
                                      sep='\t')
    
    success_rates = []
    for e, data in datasets.iteritems():
        print 'Stats for MAX_EDGE = {}'.format(e)
        exactly_1 = sum(data['count'] == 1)
        greater_than_1 = sum(data['count'] > 1)
        success_rate = (exactly_1 - greater_than_1)/len(data)
        success_rates.append(success_rate)
        print success_rate
    
    ax = yzer.set_up_plot()
    title = 'Transcript initiation recapture versus MAX_EDGE'
    yzer.add_title(title, ax)
    yzer.add_axis_labels('MAX_EDGE value', 
                         'Initiation recapture rate')
    
    yzer.plot(edges, success_rates, '-o')
    yzer.save_plot_with_dir(save_dir=img_dirpath, title=title)
    yzer.show_plot()