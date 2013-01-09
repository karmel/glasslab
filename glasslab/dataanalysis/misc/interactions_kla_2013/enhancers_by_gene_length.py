'''
Created on Jan 2, 2013

@author: karmel
'''

from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
import pandas

if __name__ == '__main__':
    enhancer_counts = True # Are we looking at enhancer interactions (False) or counts (True)?
    
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/enhancers_by_gene_length'
    dirpath = yzer.get_path(dirpath)
    
    img_dirpath = yzer.get_and_create_path(dirpath, 'scatterplots')

    counted = enhancer_counts and 'enhancer' or 'interaction'
    
    # The first set has length with interaction counts; 
    # the second has length for all transcripts, even those without interactions.
    # We want to merge such that we add the interaction-less genes with a count of 0.
    data = yzer.import_file(yzer.get_filename(dirpath,'{0}_counts_by_refseq.txt'.format(counted)))
    all_data = yzer.import_file(yzer.get_filename(dirpath,'refseq_all.txt'))
    all_data = all_data[~all_data['id'].isin(data['id'])]
    data = pandas.concat([data, all_data])
    data = data.reset_index().fillna(0)
    
    notx = data[data['sequencing_run_id'] == 765]
    kla_30m = data[data['sequencing_run_id'] == 766]
    kla_4h = data[data['sequencing_run_id'] == 773]
    no_intxns = data[data['sequencing_run_id'] == 0]
    
    # Zero won't show up in a log plot, so add one.
    no_intxns['count'] = 1
    
    
    ax = yzer.scatterplot(no_intxns, 
                     xcolname='length', ycolname='count', log=True, color='#CCCCCC', 
                     label='No {0}s'.format(counted), show_2x_range=False, plot_regression=False, 
                     show_count=False, show_correlation=False, show_legend=False, show_plot=False)
    ax = yzer.scatterplot(notx, 
                     xcolname='length', ycolname='count', log=True, color='#B5D8EB', 
                     label='Notx {0}s'.format(counted), show_2x_range=False, plot_regression=False, 
                     show_count=False, show_correlation=False, show_legend=False, show_plot=False, ax=ax)
    ax = yzer.scatterplot(kla_30m, 
                     xcolname='length', ycolname='count', log=True, color='#FFBDD8', 
                     label='KLA 30m {0}s'.format(counted), show_2x_range=False, plot_regression=False, 
                     show_count=False, show_correlation=False, show_legend=False, show_plot=False, ax=ax)
    ax = yzer.scatterplot(kla_4h, 
                     xcolname='length', ycolname='count', log=True, color='#E3AAD6', 
                     title='{0} counts as a function of gene length'.format(counted.title()), 
                     xlabel='Transcript length', ylabel='Distal {0} count'.format(counted), 
                     label='KLA 4h {0}s'.format(counted), show_2x_range=False, plot_regression=False, 
                     show_count=False, show_correlation=True, show_legend=True, 
                     save_dir=img_dirpath, show_plot=True, ax=ax)

    
    