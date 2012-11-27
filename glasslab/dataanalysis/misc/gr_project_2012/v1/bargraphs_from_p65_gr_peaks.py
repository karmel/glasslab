'''
Created on Oct 1, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from matplotlib import pyplot
import numpy

def get_filters(subdata, xcol, ycol):
    none = subdata[xcol] + subdata[ycol] == 0
    lt = (subdata[xcol] + subdata[ycol] > 0) & (subdata[xcol]*1.2 <= subdata[ycol])
    nc = (subdata[xcol] + subdata[ycol] > 0) & (subdata[xcol]*1.2 > subdata[ycol]) & (subdata[xcol] < 1.2*subdata[ycol])
    gt = (subdata[xcol] + subdata[ycol] > 0) & (subdata[xcol] >= 1.2*subdata[ycol])
    return none, lt, nc, gt

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'bargraphs_from_p65_gr')
    
    
    if True:
        for main, compare, basal_cond in (('p65','GR', 'KLA'),('GR','p65', 'Dex')):
            data = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'from_peaks', 
                                                      '{0}_kla_dex_vectors.txt'.format(main)))
        
            data = data.fillna(0)
            data = data.groupby(['id','chr_name'],as_index=False).mean()
            data = data[data['tag_count'] >= 10]
            
            
            datasets = [data[filterset] for filterset in get_filters(data, 'tag_count_3', 'tag_count_4')]
            
            width = .5
            stack_names = [s.format(compare) for s in ['No {0}', 'Less {0} in KLA+Dex', 'Has {0}, no change', 'More {0} in KLA+Dex']]
            bar_names = [s.format(main) for s in ['Less {0} in KLA+Dex', 'No change', 'More {0} in KLA+Dex']]
            xvals = numpy.arange(len(bar_names)) + width
            bottoms = numpy.array([0]*len(bar_names))
            colors = ('#E5FCC2','#9DE0AD','#45ADA8','#547980')
            totals = [len(data[filterset]) for filterset in get_filters(data, 'tag_count', 'tag_count_2')
                                            if sum(filterset) > 0]
            totals = numpy.array(totals)
            
            pyplot.figure(figsize=[10, 6])
            for i, dataset in enumerate(datasets):
                subdatasets = [dataset[filterset] for filterset in get_filters(dataset, 'tag_count', 'tag_count_2')
                                            if sum(filterset) > 0]
                yvals = numpy.array(map(len, subdatasets))/totals
                ax = pyplot.bar(xvals, yvals, width, bottom=bottoms, label=stack_names[i],color=colors[i])
                bottoms = bottoms + yvals
            
            title = 'Co-occupancy of {0} and {1} by peak change'.format(main, compare)
            pyplot.title(title)
            pyplot.ylabel('Fraction')
            pyplot.xlim([0,len(bar_names)*1.8])
            pyplot.xticks(xvals+width/2,bar_names)
            pyplot.legend()
            pyplot.savefig(yzer.get_filename(img_dirpath,title.replace(' ','_')))
            pyplot.show()