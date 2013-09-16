'''
Created on May 2, 2013

@author: karmel

Note: Made font.weight = normal and axes.titlesize = 24 in matplotlibrc
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.dataanalysis.misc.demoatlas.rpkm_to_score import PrettyAxisGrapher


if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/NAR_review_data/Post-gene'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'scatterplots')

    data = yzer.import_file(yzer.get_filename(dirpath,'post_gene_transcripts.txt'))
    refseq = yzer.import_file(yzer.get_filename(dirpath,'all_expressed_refseq.txt'))
    
    refseq_with_runoff = refseq[refseq['id'].isin(data['gene_id'])]
    refseq_no_runoff = refseq[~refseq['id'].isin(data['gene_id'])]
    if False:
        print len(refseq_no_runoff)
        print refseq_no_runoff.tail(100).to_string()
    
    # Calculate length of runoff
    data['length'] = data['transcription_end'] - data['transcription_start'] + 1
    data['gene_length'] = data['gene_end'] - data['gene_start'] + 1
    
    # What might be correlated with length of runoff?
    if False:
        yzer.scatterplot(data, 'gene_length', 'length', log=True)
        yzer.scatterplot(data, 'gene_score', 'length', log=True)
        yzer.scatterplot(data, 'score', 'length', log=True)
        yzer.scatterplot(data, 'gene_score', 'score', log=True)
        yzer.scatterplot(data, 'gene_rpkm', 'rpkm', log=True)
        yzer.scatterplot(data, 'gene_rpkm', 'length', log=True)
        yzer.scatterplot(data, 'rpkm', 'length', log=True)
        yzer.boxplot([data['gene_score'], refseq_no_runoff['score']])
        yzer.boxplot([data['gene_rpkm'], refseq_no_runoff['rpkm']])
                
        
    if True:
        for subset in (refseq_no_runoff, refseq_with_runoff):
            subset['percent_covered'] = (subset['transcription_end'] - subset['transcription_start'] + 1)\
                        /(subset['transcription_end(2)'] - subset['transcription_start(2)'] + 1)
            print 'Total: ', len(subset)
            print 'Percentage: ', len(subset)/len(refseq)
            print 'mRNA: ', sum(subset['type'] == 'mRNA'), sum(subset['type'] == 'mRNA')/len(subset)
            print 'rRNA: ', sum(subset['type'] != 'mRNA'), sum(subset['type'] != 'mRNA')/len(subset)
        # Note that some transcripts encompass two very close refseq genes; we're filtering those out for now.
        mrna_with_runoff = refseq_with_runoff[(refseq_with_runoff['type'] == 'mRNA')
                                              & (refseq_with_runoff['refseq'] == 't')
                                              & (refseq_with_runoff['percent_covered'] < 1.5)]
        mrna_no_runoff = refseq_no_runoff[(refseq_no_runoff['type'] == 'mRNA')
                                            & (refseq_no_runoff['refseq'] == 't')
                                            & (refseq_no_runoff['percent_covered'] < 1.5)]
        if False:
            yzer.histogram(mrna_with_runoff['percent_covered'])
            yzer.histogram(mrna_no_runoff['percent_covered'])
            
            
            relationships = ['is contained by', 'contains', 'overlaps with']
            with_runoff_counts = [sum(mrna_with_runoff['relationship'] == rel)
                                  for rel in relationships]
            no_runoff_counts = [sum(mrna_no_runoff['relationship'] == rel)
                                  for rel in relationships]
            yzer.piechart(with_runoff_counts, labels=relationships)
            yzer.piechart(no_runoff_counts, labels=relationships)
            

    if True:
        # Filter down to high-expression genes
        def distance_to_reg_end(row):
            if row['strand'] == 0: 
                # RefSeq annotated end - transcript end; pos if Refseq is longer 
                distance = row['transcription_end(2)'] - row['transcription_end']
            elif row['strand'] == 1: 
                # transcript start - RefSeq annotated start; pos if Refseq is longer 
                distance = row['transcription_start'] - row['transcription_start(2)']
            return distance
        for subset in (mrna_with_runoff, mrna_no_runoff):
            subset['distance'] = subset.apply(distance_to_reg_end, axis=1)
        min_score = 1
        reaches_end_with = mrna_with_runoff[(mrna_with_runoff['distance'] <= 10)
                                            & (mrna_with_runoff['score'] >= min_score)]
        reaches_end_none = mrna_no_runoff[(mrna_no_runoff['distance'] <= 10)
                                          & (mrna_no_runoff['score'] >= min_score)]
        print len(reaches_end_with), len(reaches_end_none)
        print len(reaches_end_with)+len(reaches_end_none)
        print (len(reaches_end_with)+len(reaches_end_none))/(len(mrna_with_runoff)+len(mrna_no_runoff))
        print (len(mrna_with_runoff)+len(mrna_no_runoff))
        print len(reaches_end_with)/(len(reaches_end_with)+len(reaches_end_none))
        
    if False:
        min_score = 2       
        high_exp_with = mrna_with_runoff[mrna_with_runoff['score'] >= min_score] 
        high_exp_none = mrna_no_runoff[mrna_no_runoff['score'] >= min_score] 
        
        print len(high_exp_with), len(high_exp_none)
        print len(high_exp_with)+len(high_exp_none)
        print (len(high_exp_with)+len(high_exp_none))/len(refseq)
        print len(high_exp_with)/(len(high_exp_with)+len(high_exp_none))
        
    if True:
        # Cleaned up for figures.
        
        # 5a
        yzer.scatterplot(data, 'gene_score', 'score', log=True,
                         title='Post-gene RNA Score Is Correlated\nwith the Preceding Gene Score', 
                         xlabel='Vespucci Score for RefSeq Transcript',
                         ylabel='Vespucci Score for Post-gene RNA',  
                         show_2x_range=False, plot_regression=False, square=False,
                         show_count=True, show_correlation=True, show_legend=False, 
                         save_dir=img_dirpath, show_plot=False)
        
        # S3a
        other_plot_bk = yzer.other_plot
        yzer.other_plot = PrettyAxisGrapher().other_plot
        yzer.scatterplot(data, 'gene_length', 'length', log=True,
                         title='Post-gene RNA Length Is Not Correlated\nwith the Preceding Gene Length', 
                         xlabel='RefSeq Transcript Length',
                         ylabel='Post-gene RNA Length',  
                         show_2x_range=False, plot_regression=False, set_limits=False,
                         show_count=True, show_correlation=True, show_legend=False, 
                         save_dir=img_dirpath, show_plot=False)
        yzer.other_plot = other_plot_bk
        
        # 5b
        yzer.scatterplot(data, 'gene_score', 'length', log=True,
                         title='Post-gene RNA Length Is Not Correlated\nwith the Preceding Gene Score', 
                         xlabel='Vespucci Score for RefSeq Transcript',
                         ylabel='Length of Post-gene RNA',  
                         show_2x_range=False, plot_regression=False, set_limits=False,
                         show_count=True, show_correlation=True, show_legend=False, 
                         save_dir=img_dirpath, show_plot=False)
        
        # S3b
        yzer.scatterplot(data, 'gene_rpkm', 'length', log=True,
                         title='Post-gene RNA Length Is Not Correlated\nwith the Preceding Gene RPKM', 
                         xlabel='RPKM for RefSeq Transcript',
                         ylabel='Length of Post-gene RNA',  
                         show_2x_range=False, plot_regression=False, set_limits=False,
                         show_count=True, show_correlation=True, show_legend=False, 
                         save_dir=img_dirpath, show_plot=False)
        
        # 5c
        yzer.boxplot([data['gene_score'], refseq_no_runoff['score']],
                     bar_names=['With post-gene RNA', 'Without post-gene RNA',],
                     title='Genes without Post-gene RNA\nAre Less Expressed', 
                     xlabel='RefSeq Transcript Subset',
                     ylabel='Vespucci Score',  
                     save_dir=img_dirpath, show_plot=False)
        
        # S3c
        yzer.boxplot([data['gene_rpkm'], refseq_no_runoff['rpkm']],
                     bar_names=['With post-gene RNA', 'Without post-gene RNA',],
                     title='Genes without Post-gene RNA\nAre Less Expressed (RPKM)', 
                     xlabel='RefSeq Transcript Subset',
                     ylabel='RPKM',  
                     save_dir=img_dirpath, show_plot=False)
        
    
        
        