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

def get_filters_many(subdata, xcol, ycol):
    xcol2, ycol2 = 'tag_count', 'tag_count_2'
    none = subdata[xcol] + subdata[ycol] == 0
    lt = (subdata[xcol] + subdata[ycol] > 0) & (subdata[xcol]*1.2 <= subdata[ycol]) & (subdata[xcol2] <= 1.2*subdata[ycol2])
    lt_with_gain = (subdata[xcol] + subdata[ycol] > 0) & (subdata[xcol]*1.2 <= subdata[ycol]) & (subdata[xcol2] > 1.2*subdata[ycol2])
    nc = (subdata[xcol] + subdata[ycol] > 0) & (subdata[xcol]*1.2 > subdata[ycol]) \
                & (subdata[xcol] < 1.2*subdata[ycol]) & (subdata[xcol2] <= 1.2*subdata[ycol2])
    nc_with_gain = (subdata[xcol] + subdata[ycol] > 0) & (subdata[xcol]*1.2 > subdata[ycol]) \
                & (subdata[xcol] < 1.2*subdata[ycol]) & (subdata[xcol2] > 1.2*subdata[ycol2])
    gt = (subdata[xcol] + subdata[ycol] > 0) & (subdata[xcol] >= 1.2*subdata[ycol]) & (subdata[xcol2] <= 1.2*subdata[ycol2])
    gt_with_gain = (subdata[xcol] + subdata[ycol] > 0) & (subdata[xcol] >= 1.2*subdata[ycol]) & (subdata[xcol2] > 1.2*subdata[ycol2])
    return none, lt, lt_with_gain, nc, nc_with_gain, gt, gt_with_gain

def get_filters_transcript(subdata, xcol, ycol):
    down_in_kla = subdata['kla_1_lfc'] <= -1
    nc_in_kla = subdata['kla_1_lfc'].abs() < 1
    up_in_kla = subdata['kla_1_lfc'] >= 1 & (subdata['dex_over_kla_1_lfc'] > -.58)
    trans = up_in_kla & (subdata['dex_over_kla_1_lfc'] <= -.58)
    return down_in_kla, nc_in_kla, up_in_kla, trans

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'bargraphs_from_p65_gr')
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'motifs','transcript_vectors.txt'))
    transcripts['glass_transcript_id'] = transcripts['id']
    
    if True:
        for main, compare, basal_cond in (('p65','GR', 'KLA'),('GR','p65', 'Dex')):
            data = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'from_peaks', 
                                                      '{0}_kla_dex_vectors.txt'.format(main)))
        
            data = data.merge(transcripts, how='left', on='glass_transcript_id',suffixes=['','trans'])
            data = data.fillna(0)
            #data = data.groupby(['id','chr_name'],as_index=False).mean()
            data = data[data['tag_count'] >= 10]
            
            
            datasets = [data[filterset] for filterset in get_filters(data, 'tag_count_3', 'tag_count_4')]
            
            width = .5
            stack_names = [s.format(compare, main) for s in ['No {0}', 'Less {0} in KLA+Dex', 
                                                             'Has {0}, no change', 'More {0} in KLA+Dex']]
            #stack_names = [s.format(compare, main) for s in ['No {0}', 'Less {0} in KLA+Dex', 'Less {0} in KLA+Dex and gains {1}', 
            #                                           'Has {0}, no change', 'Has {0}, no change and gains {1}', 
            #                                           'More {0} in KLA+Dex', 'More {0} in KLA+Dex and gains {1}']]
            
            bar_names = ['Down\nin KLA', 'No change\nin KLA', 'Up\nin KLA, not trans', 'Transrepressed\nby Dex']
            xvals = numpy.arange(len(bar_names)) + width
            bottoms = numpy.array([0]*len(bar_names))
            colors = ['#E5FCC2','#9DE0AD','#45ADA8','#547980','#594F4F','#5A6377','#6D5977', ]
            totals = [len(data[filterset].groupby(['id','chr_name'],as_index=False).mean()) 
                                            for filterset in get_filters_transcript(data, 'tag_count', 'tag_count_2')
                                            if sum(filterset) > 0]
            totals = numpy.array(totals)
            
            pyplot.figure(figsize=[10, 6])
            for i, dataset in enumerate(datasets):
                subdatasets = [dataset[filterset].groupby(['id','chr_name'],as_index=False).mean()
                                            for filterset in get_filters_transcript(dataset, 'tag_count', 'tag_count_2')]
                
                yvals = numpy.array(map(len, subdatasets))/totals
                ax = pyplot.bar(xvals, yvals, width, bottom=bottoms, label=stack_names[i],color=colors[i])
                bottoms = bottoms + yvals
            
            title = 'Change in {1} by transcript change'.format(main, compare)
            #title = 'Change in {0} and {1} by transcript change'.format(main, compare)
            pyplot.title(title)
            pyplot.ylabel('Fraction')
            pyplot.xlim([0,len(bar_names)*1.7])
            pyplot.ylim([0,1])
            pyplot.xticks(xvals+width/2,bar_names)
            pyplot.legend()
            pyplot.savefig(yzer.get_filename(img_dirpath,title.replace(' ','_')))
            pyplot.show()