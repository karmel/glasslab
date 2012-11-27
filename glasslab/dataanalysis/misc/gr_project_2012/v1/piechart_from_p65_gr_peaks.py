'''
Created on Oct 1, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from matplotlib import pyplot
import numpy

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'piecharts_from_p65_gr')
    
    
    if True:
        for main, compare, basal_cond in (('GR','p65', 'Dex'),('p65','GR', 'KLA'),):
            data = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'from_peaks', 
                                                      '{0}_kla_dex_vectors.txt'.format(main)))
        
            # Get nearby peaks first
            ids_with_nearby = data[(data['distance_to_tss_2'].isnull() == False) & (data['distance_to_peak_2'] <= 1000)]['id']
            
            data = data.fillna(0)
            data = data.groupby(['id','chr_name'],as_index=False).mean()
            data = data[data['tag_count'] >= 10]
        
            total = len(data)
            
            has_nearby_peak = data['id'].isin(ids_with_nearby)
            bound_by_main_not_comp_not_basal = data[~has_nearby_peak & (data['tag_count_3'] < 10)]
            bound_by_main_not_comp_basal = data[has_nearby_peak & (data['tag_count_3'] < 10)]
            bound_by_main_comp_not_basal = data[~has_nearby_peak & (data['tag_count_3'] >= 10)]
            bound_by_main_comp_basal = data[has_nearby_peak & (data['tag_count_3'] >= 10)]
            counts = [len(d) for d in (bound_by_main_not_comp_not_basal, bound_by_main_not_comp_basal,
                                       bound_by_main_comp_not_basal, bound_by_main_comp_basal)]
            labels = [s.format(main, compare, basal_cond) for s in ['Has {0} but no {1} in KLA+Dex, no {0} in {2}',
                                                                    'Has {0} but no {1} in KLA+Dex, has {0} in {2}',
                                                                    'Has {0} and {1} in KLA+Dex, no {0} in {2}',
                                                                    'Has {0} and {1} in KLA+Dex, has {0} in {2}',
                                                                    ]]
            colors = ('#E5FCC2','#9DE0AD','#45ADA8','#547980')
            
            pyplot.figure(figsize=[10,10])
            patches, texts, autotexts = pyplot.pie(counts, labels=labels, colors=colors,
                                                   autopct='%.2f%%', labeldistance=1.05)
            # Hide labels. A bit of a hack...
            [t.set_fontsize(0.01) for t in texts]
            
            title = 'Co-occupancy of {0} by {1} presence'.format(main, compare)
            pyplot.title(title)
            pyplot.legend(loc='lower left')
            pyplot.savefig(yzer.get_filename(img_dirpath,title.replace(' ','_')))
            pyplot.show()
