'''
Created on Jul 11, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'peak_scatterplots')
    
    
    if True:
        for main, compare, basal_cond in (('p65','GR', 'KLA'),('GR','p65', 'Dex')):
            data = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'from_peaks', 
                                                      '{0}_kla_dex_vectors.txt'.format(main)))
        
            data = data.fillna(0)
            data = data.groupby(['id','chr_name'],as_index=False).mean()
            
            xcolname, ycolname = 'tag_count_2', 'tag_count' #'p65_kla_tag_count', 'p65_kla_dex_tag_count',
            data = data[data[ycolname] >= 10]
            
            cond_1 = (data['tag_count_3'] == 0)
            cond_2 = (data['tag_count_3'] > 0) & (data['tag_count_3'] < data['tag_count_4'] )
            cond_3 = (data['tag_count_3'] > 0) & (data['tag_count_3'] >= data['tag_count_4'] )
            ax = None
            for show_points in (True, False):
                ax = yzer.scatterplot(data[cond_1], xcolname, ycolname,
                                    log=True, color=show_points and '#333333' or 'grey', 
                                    master_dataset=data,
                                    xlabel='{0} {1} tag count'.format(main, basal_cond), 
                                    ylabel='{0} KLA+Dex tag count'.format(main), 
                                    label='No {0} in KLA+Dex {1}'.format(compare,
                                                        show_points and ' ({0})'.format(len(data[cond_1])) or ''), 
                                    add_noise=show_points, show_points=show_points,
                                    show_2x_range=False, show_legend=False, plot_regression=(not show_points),
                                    show_count=False, show_correlation=False, set_limits=True,
                                    show_plot=False, ax=ax)
                ax = yzer.scatterplot(data[cond_2], xcolname, ycolname,
                                    log=True, color=show_points and 'blue' or 'lightblue', 
                                    master_dataset=data,
                                    xlabel='{0} {1} tag count'.format(main, basal_cond), 
                                    ylabel='{0} KLA+Dex tag count'.format(main), 
                                    label='Loses {0} in KLA+Dex {1}'.format(compare,
                                                        show_points and ' ({0})'.format(len(data[cond_2])) or ''), 
                                    add_noise=show_points, show_points=show_points,
                                    show_2x_range=False, show_legend=False, plot_regression=(not show_points),
                                    show_count=False, show_correlation=False,
                                    show_plot=False, ax=ax)
                ax = yzer.scatterplot(data[cond_3], xcolname, ycolname,
                                    log=True, color=show_points and 'red' or 'pink', 
                                    master_dataset=data,
                                    title='Peak Tag Count Comparison: {0} in {1} vs. KLA+Dex'.format(main, basal_cond),
                                    xlabel='{0} {1} tag count'.format(main, basal_cond), 
                                    ylabel='Has {0} KLA+Dex tag count'.format(main), 
                                    label='Gains/maintains {0} in KLA+Dex {1}'.format(compare, 
                                                        show_points and ' ({0})'.format(len(data[cond_3])) or ''),
                                    add_noise=show_points, show_points=show_points,
                                    show_2x_range=(not show_points), show_legend=(not show_points), 
                                    plot_regression=(not show_points),
                                    show_count=(not show_points), show_correlation=False, 
                                    show_plot=False, ax=ax)
            
                #yzer.xlim(ax, min(data[xcolname]), max(data[xcolname]))
                #yzer.ylim(ax, min(data[ycolname]), max(data[ycolname]))
            
            yzer.save_plot(yzer.get_filename(img_dirpath, '{0}_in_{1}_vs_KLA-Dex_by_{2}.png'.format(
                                                                    main, basal_cond, compare)))
            yzer.show_plot()

    if False:
        for main, compare, basal_cond in (('p65','GR', 'KLA'),('GR','p65', 'Dex')):
            data = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'from_peaks', 
                                                      '{0}_kla_dex_vectors.txt'.format(main)))
        
            data = data.fillna(1)
            data = data.groupby(['id','chr_name'],as_index=False).mean()
            
            data['ratio_main'] = data['tag_count']/data['tag_count_2']
            data['ratio_compare'] = data['tag_count_3']/data['tag_count_4']
            ax = yzer.scatterplot(data, 'ratio_main', 'ratio_compare',
                                    log=True, color='blue', 
                                    master_dataset=data,
                                    title='Ratio of {0} and {1} in basal vs. KLA+Dex'.format(main, compare),
                                    xlabel='{0} KLA+Dex/{0} {1}'.format(main, basal_cond), 
                                    ylabel='{0} KLA+Dex/{0} {1}'.format(compare, basal_cond),
                                    add_noise=False,
                                    show_2x_range=True, show_legend=False, 
                                    plot_regression=False,
                                    show_count=True, show_correlation=True, 
                                    show_plot=False)
            
            yzer.save_plot(yzer.get_filename(img_dirpath, 'ratios_comparison_{0}.png'.format(
                                                                    main, basal_cond, compare)))
            yzer.show_plot()