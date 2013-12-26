'''
Created on Dec 26, 2013

@author: karmel

From Whyte 2012:

To identify super-enhancers, we first ranked all enhancers in a cell type 
by increasing total background-subtracted ChIP-seq occupancy of Med1, 
and plotted the total background-subtracted ChIP-seq occupancy of Med1 
in units of total rpm/bp (reads per million per base pair) (Figures 1C and 4B). 
In cases where Med1 ChIP-seq data were not available, we used the total 
background subtracted ChIP-seq occupancy of the master regulator instead 
(Figure S5). These plots revealed a clear point in the distribution of 
enhancers where the occupancy signal began increasing rapidly. To geometrically 
define this point, we first scaled the data such that the x and y axis were 
from 0-1. We then found the x axis point for which a line with a slope of 1 
was tangent to the curve. We define enhancers above this point to be 
super-enhancers, and enhancers below that point to be typical enhancers.

In order to do the same here, we need to:

1. Read in merged peaks.
2. Calculate and sort by norm reads/bp.
3. Normalize reads/bp and rank.
4. Calculate the moving average for reads/bp.
5. Calculate the slope between every pair of averaged reads/bp.
6. Calculate the moving average for the slope.
7. Find the x for which the averaged slope is >= 1.
8. Output all enhancers at or above that rank.
 
'''
from __future__ import division
from pandas.io import parsers
import sys
import os
from pandas.stats.moments import rolling_mean
from matplotlib import pyplot


class SuperFinder(object):
    m = 1
    moving_window = 100
        
    def import_file(self, filename, separator='\t', header=True, 
                    index_col=None, skiprows=None):
        if header: header_row = 0
        else: header_row = None
        data = parsers.read_csv(filename, 
                                sep=separator, 
                                header=header_row, 
                                index_col=index_col,
                                skiprows=skiprows)
        
        return data
    
    def sort_by_reads_per_bp(self, data):
        '''
        Calculate reads per basepair in each peak, and sort ascending.
        Assign rank for use later.
        '''
        data['reads_per_bp'] = data['tag_count']/(data['end'] - data['start'] + 1)
        data['rank'] = data['reads_per_bp'].rank(method='first') 
        data = data.sort('rank')
        return data
    
    def normalize(self, data):
        '''
        Normalize both reads_per_bp and rank such that the range is 0 to 1.
        '''
        for col in ('reads_per_bp', 'rank'):
            data[col + '_norm'] = (data[col] - data[col].min())/ data[col].max()
            
        return data
    
    
    def moving_average(self, data, column='reads_per_bp_norm'):
        '''
        Calculate moving average of slope over defined window.
        '''      
        data['avg_' + column] = rolling_mean(data[column], 
                                            self.moving_window)
        
        return data
    
    def pairwise_slopes(self, data):
        '''
        Calculate the slope between every pair of points (rank, reads_per_bp)
        '''
        data['slope'] = 0.0
        for x in xrange(len(data) - 1):
            # Convert the row-number to row-index so that we can select
            # from the dataframe the rows we want in order.
            idx, next_idx = data.index[x], data.index[x+1]
            run = data.ix[next_idx]['rank_norm'] - data.ix[idx]['rank_norm']
            rise = data.ix[next_idx]['avg_reads_per_bp_norm'] \
                    - data.ix[idx]['avg_reads_per_bp_norm']
            data.set_value(idx, 'slope', rise/run)
        
        return data
    

    def find_breakpoint(self, data):
        '''
        Given the moving average, find the point at which it passes
        the slope of interest. Return that rank. 
        '''
        # We want to take the last rank for which the slope passes m
        slope_small = data[data['avg_slope'] < self.m]
        min_rank = slope_small['rank_norm'].max()
        supers = data[data['rank_norm'] >= min_rank]
        return supers
    
    def output_supers(self, data, filename):
        output_file, _ = os.path.splitext(filename)
        output_file = output_file + '_supers.txt'
        
        data = data[['chr', 'start', 'end', 
                     'reads_per_bp', 'tag_count', 
                     'merged_ids']]
        data = data.sort(['reads_per_bp'], ascending=False)
        
        data.to_csv(output_file, sep='\t',  
                    header=True, index=True)


    def plot(self, data, supers, filename):
        output_file, _ = os.path.splitext(filename)
        output_file = output_file + '_supers.png'
        
        pyplot.plot(data['rank_norm'], data['reads_per_bp_norm'], 
                    label='Raw')
        pyplot.plot(data['rank_norm'], data['avg_reads_per_bp_norm'], 
                    label='Moving Avg (window = {})'.format(self.moving_window))
        
        min_rank = supers['rank_norm'].min()
        pyplot.plot([min_rank, min_rank], [0,1], label='Cut-off')
        
        pyplot.title('Super-enhancer Indentification')
        pyplot.xlabel('Rank (norm)')
        pyplot.ylabel('Reads per million per basepair (norm)')
        pyplot.legend()
        
        pyplot.savefig(output_file)
        pyplot.show()
        
    def find_superenhancers(self, filename):
        data = finder.import_file(filename, index_col=0)
        data = finder.sort_by_reads_per_bp(data)
        data = finder.normalize(data)
        data = finder.moving_average(data, 'reads_per_bp_norm')
        data = finder.pairwise_slopes(data)
        data = finder.moving_average(data, 'slope')
        supers = finder.find_breakpoint(data)
        finder.output_supers(supers, filename)
        finder.plot(data, supers, filename)
    
if __name__ == '__main__':
    
    # Get filepath
    try: filename = sys.argv[1]
    except IndexError: filename = '/Users/karmel/GlassLab/Notes_and_Reports/Super-Enhancers/whyte_2012/th1_tbet_peaks_merged.txt'
    
    finder = SuperFinder()
    finder.find_superenhancers(filename)
    
