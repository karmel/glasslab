'''
Created on Oct 24, 2013

@author: karmel

We want to evaluate the Hah et al HMM using the GRO5 data from 
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
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/NAR_review_data_2/hmm_genome_wide'
    dirpath = yzer.get_path(dirpath)
    #img_dirpath = yzer.get_and_create_path(dirpath, 'success_rates')
    
    datasets = OrderedDict()
    for suffix in ('dm_10000','hmm_100_5'):
        datasets[suffix] = pandas.read_csv(yzer.get_filename(dirpath, 
                                    'gro5_count_{}.txt'.format(suffix)),
                                      header=0, 
                                      sep='\t')
    
    success_rates = []
    for suffix, data in datasets.iteritems():
        print 'Stats for {}'.format(suffix)
        exactly_1 = sum(data['count'] == 1)
        greater_than_1 = sum(data['count'] > 1)
        success_rate = (exactly_1 - greater_than_1)/len(data)
        success_rates.append(success_rate)
        print success_rate
    
    