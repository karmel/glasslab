'''
Created on Mar 20, 2012

@author: karmel

Miscellaneous methods for graphing tag count plots.
'''
from __future__ import division

import math
from string import capwords

from matplotlib import pyplot, font_manager
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import ScalarFormatter
import numpy
from numpy.lib.polynomial import poly1d

from glasslab.dataanalysis.base.datatypes import TranscriptAnalyzer
class SeqGrapher(TranscriptAnalyzer):
    fig_size = 10
    disable_show_plot = False
    legend_location = 'upper left'
    legend_columns = 1
    
    def scatterplot(self, data, xcolname, ycolname, 
                    subplot=111, log=False, color='blue',
                    master_dataset=None,
                    title='', xlabel=None, ylabel=None,
                    label='', add_noise=False,
                    show_points=True, set_limits=False, square=False,
                    show_2x_range=True, plot_regression=False,
                    show_count = True, show_correlation=True,
                    text_shift=0, text_color='black', 
                    show_legend=True, save_dir='', save_name='', 
                    show_plot=True, ax=None):
        
        '''
        Designed to show scatterplots of tags by tags for a given run
        
        Sample using two colors: 
        # Split into those that are also diff in non-diabetic NOD and those not
        refseq_diabetic = refseq[abs(refseq['balb_nod_notx_1h_fc']) < 1]
        refseq_strain = refseq[abs(refseq['balb_nod_notx_1h_fc']) >= 1]
        
        ax = grapher.scatterplot(refseq_diabetic, xcolname, ycolname,
                            log=True, color='blue', master_dataset=refseq,
                            title='Nonplated Diabetic NOD vs. BALBc Refseq Transcripts',
                            xlabel='BALBc notx', ylabel='NOD notx', label='Different only with diabetes when not plated',
                            show_2x_range=False, show_legend=False, 
                            show_count=False, show_correlation=False, show_plot=False)
        ax = grapher.scatterplot(refseq_strain, xcolname, ycolname,
                            log=True, color='red', master_dataset=refseq,
                            label='Different in NOD without diabetes (plated)',
                            show_2x_range=True, show_legend=True,
                            show_count=True, show_correlation=True, show_plot=False)
        grapher.save_plot(os.path.join(dirpath, 'nonplated_diabetic_nod_vs_balbc_scatterplot.png'))
        grapher.show_plot()
        '''
        # Sometimes we want to use a superset of the data for plot setup.
        master_dataset = master_dataset or data
        
        ax = self.set_up_plot(ax, subplot)
        
        if set_limits:
            self.xlim(ax, min(data[xcolname]), max(data[xcolname]))
            self.ylim(ax, min(data[ycolname]), max(data[ycolname]))
        if square:
            min_val = min(min(data[xcolname]), min(data[ycolname]))
            max_val = max(max(data[xcolname]), max(data[ycolname]))
            self.xlim(ax, min_val, max_val)
            self.ylim(ax, min_val, max_val)
        # Add some noise to prevent overlap?
        if add_noise:
            xcol = data[xcolname]*(1 + .01*numpy.random.randn(len(data[xcolname])))
            ycol = data[ycolname]*(1 + .01*numpy.random.randn(len(data[ycolname])))
        else: xcol, ycol = data[xcolname], data[ycolname]
        
        if show_points:
            # Plot points
            pyplot.plot(xcol, ycol,
                        'o', markerfacecolor='None',
                        markeredgecolor=color, label=label)
        
        # Log scaled?
        if log:
            if log <= 1: log = 2
            pyplot.xscale('log', basex=log, nonposx='clip')
            pyplot.yscale('log', basey=log, nonposy='clip')
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_major_formatter(ScalarFormatter())
        
        # Labels
        if xlabel is None: xlabel = capwords(xcolname.replace('_',' '))
        if ylabel is None: ylabel = capwords(ycolname.replace('_',' '))
                
        self.add_axis_labels(xlabel, ylabel)
        
        self.add_title(title, ax)
        
        # Plot a least squares linear regression?
        if plot_regression:
            x, y = data[xcolname].copy(), data[ycolname].copy()
            # Reduce Pandas Series to Numpy arrays to avoid indexing problems, if necessary
            try: x, y = x.values, y.values
            except AttributeError: pass
            y = y[x.argsort()]
            x.sort()
            m, c = numpy.polyfit(x, y, 1)
            pyplot.plot(x, m*x+c, color=color, label=label + ' (Linear fit)')
        
        # Show lines at two-fold change?
        if show_2x_range:
            if show_2x_range is True:  fold_change= 2.0
            else: fold_change = float(show_2x_range)
            pyplot.plot([1, (1/fold_change)*max(master_dataset[ycolname])], 
                        [fold_change, max(master_dataset[ycolname])], 
                        '--', color='black', label='{0}-fold change'.format(fold_change))
            pyplot.plot([1, max(master_dataset[xcolname])], [(1/fold_change), (1/fold_change)*max(master_dataset[xcolname])], 
                        '--', color='black')
        
        # Show Pearson correlation and count?
        if text_shift is True: text_shift = .1
        elif text_shift is False: text_shift = 0
        if text_color is True: text_color=color
        elif text_color is False: text_color = 'black'
        if show_count: self.show_count_scatterplot(master_dataset, ax, text_shift, text_color)
        if show_correlation: self.show_correlation_scatterplot(master_dataset, xcolname, ycolname, ax, text_shift, text_color)
        
        if show_legend: self.legend()
            
            
        # Any other operations to tack on?
        self.other_plot(ax=ax)
        
        self.save_plot_with_dir(save_dir, save_name, title)
        if show_plot: self.show_plot()
    
        return ax
    
    def set_up_plot(self, ax=None, subplot=111, wide=False):
        # Set up plot
        width = self.fig_size
        if wide: width = width*(wide is True and 2 or wide)
             
        if not ax: pyplot.figure(figsize=[width*int(str(subplot)[1]), # Width for cols
                                          self.fig_size*int(str(subplot)[0])]) # Height for rows
        ax = pyplot.subplot(subplot)
        return ax
        
    def show_count_scatterplot(self, data, ax, text_shift=0, text_color='black'):
        # Scale for font size-- move left to avoid spilling over right axis
        size_scalar = min(1, 24/font_manager.FontManager().get_default_size())
        pyplot.text(size_scalar*.75, .125+text_shift, 'Total count: %d' % len(data),
                        color=text_color,
                        transform=ax.transAxes)
        
    def show_correlation_scatterplot(self, data, xcolname, ycolname, ax, text_shift=0, text_color='black'):
        # Scale for font size-- move left and down to avoid overlap with count
        size_scalar = min(1, 24/font_manager.FontManager().get_default_size())
        pyplot.text(size_scalar*.75, size_scalar*(.1+text_shift), 'r = %.4f' % data[xcolname].corr(data[ycolname]),
                    color=text_color,
                    transform=ax.transAxes)
            
    def show_plot(self): 
        if not self.disable_show_plot: 
            pyplot.show()
        pyplot.close()
    def save_plot(self, filename): pyplot.savefig(filename)
    def save_plot_with_dir(self, save_dir='', save_name='', title=''):
        if save_dir: 
            self.save_plot(self.get_filename(save_dir, 
                        save_name or (title.replace(' ','_').replace('\n','_') + '.png')))
        
    
    def add_title(self, title, ax):
        pyplot.title(title)
        ax.title.set_y(1.02)   
        
    def add_axis_labels(self, xlabel, ylabel):
        pyplot.xlabel(xlabel or '', labelpad=10)
        pyplot.ylabel(ylabel or '', labelpad=10)
        
    def xlim(self, ax, min_x, max_x):
        ax.set_xlim([min_x, max_x])

    def ylim(self, ax, min_y, max_y):
        ax.set_ylim([min_y, max_y])
    
    def legend(self):
        pyplot.legend(loc=self.legend_location, ncol=self.legend_columns)
        
    def other_plot(self, *args, **kwargs):
        return True
    
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
            ascent = [x/(segments+1) for x in range(0,segments+2)][1:-1]
            colors_r = [1]*segments + ascent[::-1] + [0]*segments
            colors_g = ascent + [1]*segments + ascent[::-1]
            colors_b = [0]*segments + ascent + [1]*segments
            base_cols = zip(colors_r, colors_g, colors_b)
            
        # Take from beginning and end of list preferentially
        selected = base_cols[:max(1,int(math.ceil(number/2)))] + base_cols[-int(number/2):]
        return selected

    
    def bargraph_for_transcript(self, transcript_row, cols,
                                bar_names=None, 
                                title='', xlabel='', ylabel='',
                                convert_log_to_fc=True, show_2x_range=True,
                                show_plot=True):
        '''
        Designed to draw bar graphs for fold changes across many
        runs for a single transcript.
        
        Sample for creation of bar graph:
        gene_row = refseq[refseq['gene_names'] == '{Clec4e}']
        grapher.bargraph_for_transcript(gene_row, 
                                        ['balb_nod_notx_1h_fc', 'balb_nod_kla_1h_fc',
                                         'diabetic_balb_nod_notx_1h_fc', 'diabetic_balb_nod_kla_1h_fc',
                                         'nonplated_diabetic_balb_nod_notx_fc',],
                                        bar_names=['Non-diabetic\nnotx 1h', 'Non-diabetic\nKLA 1h',
                                                   'Diabetic\nnotx 1h', 'Diabetic\nKLA 1h',
                                                   'Nonplated diabetic\nnotx 1h',],
                                        title='Mincle (Clec4e) Fold Change in NOD vs. BALBc GRO-seq',
                                        ylabel='Fold Change in NOD vs. BALBc',
                                        show_plot=False)
        grapher.save_plot(os.path.join(dirpath, 'clec4e_fold_change_bargraph.png'))
        grapher.show_plot()
        '''
        
        axis_spacer = .2
        bar_width = .8
        xvals = [x + axis_spacer for x in range(0,len(cols))]
        
        default_ylabel = 'Fold Change'
        if convert_log_to_fc:
            yvals = [2**transcript_row[col] for col in cols]
        else:
            yvals = [transcript_row[col] for col in cols]
            default_ylabel = 'Log(2) ' + default_ylabel
        
        # Set up plot
        ax = pyplot.subplot(111)
        ax.set_xlim([0,len(cols) + axis_spacer])
        # Show line at zero
        pyplot.plot([0,len(cols) + axis_spacer], [0,0], '-',color='black')
        
        # And at two-fold change
        if show_2x_range and max(yvals) >= 2:
            pyplot.plot([0,len(cols) + axis_spacer], 
                        [convert_log_to_fc and 2 or 1, convert_log_to_fc and 2 or 1, ], 
                        '--',color='black')
        if show_2x_range:
            pyplot.plot([0,len(cols) + axis_spacer], 
                        [convert_log_to_fc and 1 or 0, convert_log_to_fc and 1 or 0, ], 
                        '-',color='#CCCCCC')
            pyplot.plot([0,len(cols) + axis_spacer], 
                        [convert_log_to_fc and .5 or -1, convert_log_to_fc and .5 or -1, ], 
                        '--',color='black')
            
        ax.bar(xvals, yvals, bar_width, color='#C9D9FB')
        
        ax.set_xticks([x + .5*bar_width for x in xvals])
        ax.set_xticklabels(bar_names or cols)
        
        ylabel = ylabel or default_ylabel
        self.add_axis_labels(xlabel or 'Experimental Conditions', ylabel)
        self.add_title(title, ax)
        
        # Any other operations to tack on?
        self.other_plot()
        
        if show_plot: self.show_plot()
    
        return ax
    
    def boxplot(self, data,
                    bar_names=None, subplot=111, wide=False,
                    title='', xlabel=None, ylabel=None,
                    show_outliers=False, save_dir='', save_name='', 
                    show_plot=True, ax=None):
        '''
        Draw a boxplot for passed data.
        
        Data should be a list of lists of vals same length as bar_names::
            [[1,2,3,4,1,2,3,4],
            [4,5,6,7,4,5,6,7],
            ...]
        '''
        ax = self.set_up_plot(ax, subplot, wide=wide)
        
        if show_outliers: 
            symbol = '+'
        else: 
            symbol = ''
        bp = pyplot.boxplot(data, sym=symbol)
        
        if bar_names: 
            ax.set_xticklabels(bar_names)
        
        self.add_axis_labels(xlabel, ylabel)
        self.add_title(title, ax)
        #pyplot.ylim(0,100)
        # Any other operations to tack on?
        self.other_plot()
        
        self.save_plot_with_dir(save_dir, save_name, title)
        if show_plot: self.show_plot()
        
        return ax
    
    def histogram(self, data,
                    subplot=111, bins=50,
                    title='', xlabel=None, ylabel=None, show_legend=True,
                    label='', fill=True,
                    color=None, save_dir='', save_name='',
                    show_plot=True, ax=None):
        '''
        Draw a histogram for passed data.
        
        Data should be a list of vals::
            [1,2,3,4,1,2,3,4]
        
        '''
        ax = self.set_up_plot(ax, subplot)
        
        bp = pyplot.hist(data, bins=bins, color=color or self.get_colors(1)[0], 
                         fill=fill, label=label)
        
        self.add_axis_labels(xlabel, ylabel)
        self.add_title(title, ax)
        if show_legend: self.legend()
        
        #pyplot.ylim(0,100)
        # Any other operations to tack on?
        self.other_plot()
        
        self.save_plot_with_dir(save_dir, save_name, title)
        if show_plot: self.show_plot()
        
        return ax
    
    def timeseries(self, dates, data, subplot=111,
                   show_average=False, show_median=True, 
                   colors=None, labels=None,
                   title='', xlabel=None, ylabel=None,
                   show_legend=True, show_plot=True, ax=None):
        
        '''
        Data expected as a list of lists::
        
            [[group 1 member 1 val 1, group 1 member 1 val 2, group 1 member 1 val 3 ...],
                [group 1 member 2 val 1, group 2 member 1 val 2, group 1 member 2 val 3 ...]],
            [[group 2 member 1 val 1, group 2 member 1 val 2, group 2 member 1 val 3 ...] ...
            
        '''
        ax = self.set_up_plot(ax, subplot)
        
        
        for i, group in enumerate(data):
            color = colors and colors[i] or None
            
            for j, member in enumerate(group):
                # Only include label with first member
                if j == 0: label = labels and labels[i] or None
                else: label = None
            
                ax.plot(dates, member, 'o', markerfacecolor='None',
                        markeredgecolor=color, label=label)
        
        # Go through loops again to group by type (point, avg, median) in legend 
        if show_average:
            for i, group in enumerate(data):
                avg = numpy.array(group).mean(axis=0)
                ax.plot(dates, avg, '-', color=colors[i], label='Avg {0}'.format(labels[i]))
     
        store_median = None
        if show_median:
            for i, group in enumerate(data):
                med = numpy.median(group, axis=0)
                if store_median is None: store_median = med
                ax.plot(dates, med, '-', color=colors[i], label='Median {0}'.format(labels[i]))
        
        ax.plot(dates, (store_median - med) + 100, '-', color='green', linewidth=3, label="Karmel's level of excitement")
    
        from datetime import datetime, timedelta
        ax.plot([datetime(2012,4,3)]*2, [0,600], '--', color='black', label='IP injection')
        ax.plot([dates[0] - timedelta(days=1),dates[-1:][0] + timedelta(days=1)], 
                [200,200], '-', color='#CCCCCC', label='Threshold for diabetes')
        
        ax.figure.autofmt_xdate()
        self.add_axis_labels(xlabel, ylabel)
        self.add_title(title, ax)
        
        if show_legend: self.legend()
        
        # Any other operations to tack on?
        self.other_plot()
        
        if show_plot: self.show_plot()
        
        return ax
    
    def piechart(self, counts, labels,
                 title='', colors=None, small_legend=False,
                 save_dir='', save_name='', show_plot=True):
        
        ax = self.set_up_plot()
        spectrum = self.get_colors(len(counts))
        patches, texts, autotexts = pyplot.pie(counts, labels=labels, 
                                               colors=colors or spectrum,
                                               autopct='%.2f%%')
        # Hide labels. A bit of a hack...
        [t.set_text('') for t in texts]
        
        pyplot.title(title + ' ({0} total)'.format(sum(counts)))
        
        font_p = FontProperties()
        if small_legend:
            font_p.set_size(10)
        self.legend()
        
        self.save_plot_with_dir(save_dir, save_name, title)
        if show_plot: self.show_plot()
    
    def bargraph_for_transcripts(self, data, indices=None, cols=None,
                                 bar_names=None, 
                                 title='', xlabel='', ylabel='', rank_label='',
                                 convert_log_to_fc=True, show_2x_range=True,
                                 show_plot=True, save_dir='', save_name='', ):
        '''
        For multiple transcripts, show bar graphs, including rank marker.
        '''
        bar_count = len(indices or data)
        axis_spacer = .4
        bar_width = max(8/bar_count, 1.6)
        xvals_per_col = [x*1.25*bar_width + axis_spacer for x in range(0,bar_count)]
        
        # Get bar positions for all columns, adding an extra space between groups 
        xvals = [map(lambda x: x*col_number + axis_spacer*(col_number - 1), xvals_per_col)
                    for col_number in range(1,len(cols or [1])+1)] 
        # Unnest
        xvals = [item for sublist in xvals for item in sublist]
        
        fig = pyplot.figure(figsize=[max(len(xvals)*bar_width,self.fig_size), self.fig_size])
        # Set up plot
        ax = pyplot.subplot(111)
        max_x = max(xvals) + axis_spacer + bar_width
        # Show line at zero
        pyplot.plot([0,max_x], [0,0], '-',color='black')
        default_ylabel = 'Fold Change'
        
        # And at two-fold change
        if show_2x_range:
            pyplot.plot([0,max_x], 
                        [convert_log_to_fc and 2 or 1, convert_log_to_fc and 2 or 1, ], 
                        '--',color='black')
            pyplot.plot([0,max_x], 
                        [convert_log_to_fc and 1 or 0, convert_log_to_fc and 1 or 0, ], 
                        '-',color='#CCCCCC')
            pyplot.plot([0,max_x], 
                        [convert_log_to_fc and .5 or -1, convert_log_to_fc and .5 or -1, ], 
                        '--',color='black')
            
        if cols:
            # Whole dataset passed; process
            yvals = []
            for col in cols:
                if convert_log_to_fc:
                    yvals += [2**data.ix[index][col] for index in indices]
                else:
                    yvals += [data.ix[index][col] for index in indices]
                    default_ylabel = 'Log(2) ' + default_ylabel
        else: yvals = data
                        
        ax.bar(xvals, yvals, bar_width, color='#C9D9FB')
        
        pyplot.ylabel(ylabel or default_ylabel)
        #pyplot.xlabel(xlabel or 'Gene')
        
        bar_middles = [x + .5*bar_width for x in xvals]
        # Add ranks along right axis
        if indices:
            try: 
                r_ax = ax.twinx()
                r_ax.plot(bar_middles, data.ix[indices]['rank'],'v',markerfacecolor='None',markeredgecolor='blue')
                r_ax.set_ylabel(rank_label or 'Rank')
                
            except KeyError: pass
            
        ax.set_xticks(bar_middles)
        ax.set_xticklabels(bar_names or cols)
        ax.set_xlim([0,max_x])
        
        
        
        self.add_title(title, ax)
        
        # Any other operations to tack on?
        self.other_plot()
        
        self.save_plot_with_dir(save_dir, save_name, title)
        if show_plot: self.show_plot()
    
        return ax
    
    def plot(self, *args, **kwargs):
        '''
        Just pass to pyplot.
        '''
        return pyplot.plot(*args, **kwargs)