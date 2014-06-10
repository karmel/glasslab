'''
Created on Dec 27, 2013

@author: karmel

Create chromosome-by-chromosome ideograms with peaks plotted from each passed
file so that we can compare superenhancer calls with different peak types.

'''
from __future__ import division
import sys
import os
from matplotlib import pyplot
from glasslab.superenhancers.base import PeakFinder
from matplotlib.collections import BrokenBarHCollection
import math


class ChromosomePlotter(PeakFinder):
    datasets = None
    height = 1
    spacing = 0.2
    names = None
    path = ''
    xranges = None
    
    def __init__(self):
        self.datasets = []
        self.names = []
        self.xranges = {}
        
    def import_files(self, filenames):
        for f in filenames:
            data = self.import_file(f, header=True, index_col=False)
            self.datasets.append(data)
    
    def set_names(self, filenames):
        '''
        Determine peak names from filenames; assumes that peak name is before
        the first underscore in the passed file basename.
        
        The last filename bequeaths the path for saving.
        '''    
        for f in filenames:
            basename = os.path.basename(f)
            peak_name = basename.split('_')[0]
            self.names.append(peak_name)
            self.path = os.path.dirname(f)
        
    def determine_widths(self):
        for i, data in enumerate(self.datasets):
            self.datasets[i]['width'] = data['end'] - data['start'] + 1
            
    def chromsome_xrange(self, chr_data):
        '''
        Get xranges (positions) for passed peaks for a single chromosome.
        '''
        return chr_data[['start','width']].to_records(index=False).tolist()
    
    def xranges_for_each_chr(self, data):
        '''
        For each chromosome in the passed data, get xrange for plotting.
        '''
        grouped = data.groupby('chr')
        for name, group in grouped:
            try: 
                self.xranges[name].append(self.chromsome_xrange(group))
            except KeyError:
                self.xranges[name] = [self.chromsome_xrange(group)]
                
    def determine_xranges(self):
        for data in self.datasets:
            self.xranges_for_each_chr(data)
        
    def plot_bars(self, ax):
        '''
        Once xranges are determined, create bars for each chr and sample.
        '''
        
        # Some fancy sorting to get the chromosomes in 1, 2. ..10, M order
        keys = sorted(self.xranges, 
                      key=lambda x:x.replace('chr','').isdigit() \
                        and x.replace('chr','').rjust(2,'0') \
                        or x.replace('chr',''))
        
        ystart = 0
        colors = self.get_colors(len(self.datasets))
        yticks = []
        yticklabels = []

        for chrom in keys:
            for i, xranges in enumerate(self.xranges[chrom]):
                yrange = (ystart, self.height)
                
                coll = BrokenBarHCollection(xranges, yrange, color=colors[i])
                ax.add_collection(coll)
                
                # save ticks and labels for plotting
                center = yrange[0] + yrange[1]/2.
                yticks.append(center)
                
                # Only label the topmost row with the chromosome.
                label = i == (len(self.xranges[chrom]) - 1) and chrom or ''
                yticklabels.append(label)
                
                # Increment ystart for next iteration
                ystart += self.height + self.spacing 

        return yticks, yticklabels

    def get_colors(self, number):
        '''
        Set up colors so that we can set a range from red to green
        Hold R at 1 while we pass through yellow, then decrement, then hold at 0 through blue;
        meanwhile, increment G while we pass through yellow, then hold at 1 through greens, then decrement to blue;
        and, hold B at 0 through yellow and green, then increment til blue.
        
        Add grey to the front for background.
        '''
        if number < 7: 
            base_cols = ('#FFABAB','#FFDAAB','#FFFB97','#DDFFAB','#ABE4FF','#D9ABFF')
        else:
            # Split needed number into three; for first seg, hold red steady.
            # Then green for the second segment, then blue.
            segments = max(int(math.ceil(number/3)),1)
            ascent = [x/(segments+1) for x in xrange(0,segments+2)][1:-1]
            colors_r = [1]*segments + ascent[::-1] + [0]*segments
            colors_g = ascent + [1]*segments + ascent[::-1]
            colors_b = [0]*segments + ascent + [1]*segments
            base_cols = zip(colors_r, colors_g, colors_b)
            
        # Take from beginning and end of list preferentially
        selected = base_cols[:max(1,int(math.ceil(number/2)))] + base_cols[-int(number/2):]
        return selected

    def plot_peaks(self):
        fig = pyplot.figure(figsize=[10,10*len(self.datasets)])
        ax = fig.add_subplot(111)
        
        yticks, yticklabels = self.plot_bars(ax)
        ax.axis('tight')
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.set_xticks([])
        
        # Add key at the bottom of the plot
        colors = self.get_colors(len(self.names))
        offset = .02
        for i, name in enumerate(self.names):
            pyplot.figtext(offset + .1, .08, name, 
                           color=colors[i], weight='bold')
            offset += .2

        ax.set_title('Chromosomal Locations of Super-enhancers')
        save_name = '_'.join(self.names) + '_chr_locations.png'
        pyplot.savefig(os.path.join(self.path, save_name))
        pyplot.show()
        
    def draw_plots(self, filenames):
        '''
        Put it all together.
        '''
        self.import_files(filenames)
        self.set_names(filenames)
        self.determine_widths()
        self.determine_xranges()
        self.plot_peaks()
        
if __name__ == '__main__':
    
    # Get filepath
    if len(sys.argv) > 1: filenames = sys.argv[1:]
    else:
        filenames = ['/Users/karmel/GlassLab/Notes_and_Reports/Super-Enhancers/whyte_2012/th1_tbet_peaks_merged_tag_count_supers.txt',
                     '/Users/karmel/GlassLab/Notes_and_Reports/Super-Enhancers/whyte_2012/th1_tbet_peaks_merged_tag_count_supers.txt',
                     ]
    plotter = ChromosomePlotter()
    plotter.draw_plots(filenames)
    