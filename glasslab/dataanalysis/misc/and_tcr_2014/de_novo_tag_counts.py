'''
Created on Jun 8, 2014

@author: karmel

Do de novo enhancers in PCC have more tags than those in K99A or NoPep?
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/AND_TCR/Analysis/Chips1_2/Motifs'
    dirpath = yzer.get_path(dirpath)
    
    savepath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/AND_TCR/Analysis/Chips1_2/Figures/De_novo_boxplots'
    savepath = yzer.get_path(savepath)
    
    for ab in ('me2', 'ac'):
        conditions = ('K99A','PCC')
        all_tag_counts = []
        de_novo_tag_counts = []
        for peptide in conditions:
            pep_dirpath = yzer.get_filename(dirpath, 
                            '{}_{}'.format(peptide, ab))
            
            
            filename = yzer.get_filename(pep_dirpath, 
                            '{}_{}_enhancers_batf.txt'.format(peptide, ab))
            
            data = yzer.import_file(filename)
            data = data.fillna(0)
            subset = data[data['no_pep_tag_count'] == 0]
            #subset = data[data['tag_count(2)'] == 0]
            #subset = subset[subset['tag_count(3)'] == 0]
            all_tag_counts.append(data['tag_count'])
            de_novo_tag_counts.append(subset['tag_count'])
        
        # Plot as boxplot
        ax = yzer.boxplot(all_tag_counts, conditions, 
                          title='Tag Counts in {} Enhancers'.format(ab.title()), 
                          xlabel='Condition', 
                          ylabel='Normalized tag count',
                          save_dir=savepath, 
                          show_plot=True)
        ax = yzer.boxplot(de_novo_tag_counts, conditions, 
                          title='Tag Counts in de novo {} Enhancers'.format(ab.title()), 
                          xlabel='Condition', 
                          ylabel='Normalized tag count',
                          save_dir=savepath, 
                          show_plot=True)
            