'''
Created on Oct 1, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_from_p65_gr')
    
    
    if True:
        for main, compare, basal_cond in (('p65','GR', 'KLA'),('GR','p65', 'Dex')):
            data = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'from_peaks', 
                                                      '{0}_kla_dex_vectors.txt'.format(main)))
        
            data = data.fillna(0)
            data = data.groupby(['id','chr_name'],as_index=False).mean()
            data = data[data['tag_count_2'] > 0]
            
            colname = 'tag_count_diff'
            data[colname] = (data['tag_count'] - data['tag_count_2'])/data['tag_count']
            
            cond_1 = (data['tag_count_3'] == 0)
            cond_2 = (data['tag_count_3'] > 0) & (data['tag_count_3'] < data['tag_count_4'] )
            cond_3 = (data['tag_count_3'] > 0) & (data['tag_count_3'] >= data['tag_count_4'] )
            
            title = 'Difference in {0} peak tag counts by {1}'.format(main, compare)
            names = [s.format(compare) for s in ['No {0} in KLA+Dex','Loses {0} in KLA+Dex','Gains/maintains {0} in KLA+Dex']]
            ax = yzer.boxplot([data[cond_1][colname], data[cond_2][colname], data[cond_3][colname]], 
                         names,
                         title=title, 
                         xlabel='Condition', 
                         ylabel='{0} KLA+Dex tags in peak - {0} {1} tags in peak'.format(main, basal_cond), 
                         show_outliers=False, show_plot=False)
            yzer.save_plot(yzer.get_filename(img_dirpath, title.replace(' ','_')))
            yzer.show_plot()