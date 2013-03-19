'''
Created on Mar 7, 2013

@author: karmel

Quick script to plot diabetes diagnosis in different treatment groups.
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero
from matplotlib import pyplot

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/NOD_BALBc/ThioMacs/Analysis_2013_02/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'diabetes_rates')
    
    # Start from week 5; then skip to week 11
    x_labels = [5] + range(11,20)
    x_vals = range(len(x_labels))
    ctl_y = [100, 100, 100, 100, 87.5, 87.5, 87.5, 87.5, 75, 75, 75, 75, 62.5]
    ctl_x = [0, 1, 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 9]
    low_y = [100, 100, 75, 75, 75, 75, 75, 50, 50, 50, 50, 50]
    low_x = [0, 1, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9]
    med_y = [100, 100, 100, 75, 75, 50, 50, 50, 50, 25, 25, 25, 25]
    med_x = [0, 1, 2, 2, 3, 3, 4, 5, 6, 6, 7, 8, 9]
    high_y = [100, 100, 100, 100, 75, 75, 50, 50, 50, 50, 50, 50, 25]
    high_x = [0, 1, 2, 3, 3, 4, 4, 5, 6, 7, 8, 9, 9]
    
    ax = yzer.set_up_plot()
    yzer.add_axis_labels('Week', 'Percent Normoglycemic')
    title = 'Diabetes Induction with In Vivo TDB Treatment'
    yzer.add_title(title, ax)
    
    pyplot.plot(ctl_x, ctl_y, 'o-', color='black', label='Control', linewidth=5)
    pyplot.plot(low_x, low_y, 'o-', color='blue', label='Low Dose TDB (107 ug)', linewidth=4)
    pyplot.plot(med_x, med_y, 'o--', color='green', label='Medium Dose TDB (214 ug)', linewidth=4)
    pyplot.plot(high_x, high_y, 'o--', color='red', label='High Dose TDB (321 ug)', linewidth=3)
    
    pyplot.legend(loc='bottom left')
    pyplot.xticks(x_vals, x_labels)
    yzer.ylim(ax, 0, 100)
    yzer.save_plot_with_dir(img_dirpath, None, title)
    yzer.show_plot()
    
    
    