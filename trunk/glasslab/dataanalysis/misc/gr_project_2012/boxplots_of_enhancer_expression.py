'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from matplotlib import pyplot

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_by_expression')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_gene.txt'))
    data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
    
    draw_pies = True
    min_tags = 30
    ratio = 1.5
    # Make sure we have dimethyl
    data = data[data.filter(like='h3k4me2').max(axis=1) > min_tags]
    data = data[data['minimal_distance'] >= 1000]
    
    #data = yzer.collapse_strands(data)
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    transcripts['nearest_refseq_transcript_id'] = transcripts['id']
    data = data.merge(transcripts, how='left', on='nearest_refseq_transcript_id', suffixes=['','_trans'])
    
    data = data.fillna(0)
    
    transrepressed = data[(data['kla_1_lfc_trans'] >= 1) & (data['dex_over_kla_1_lfc_trans'] <= -.58)]
    not_trans = data[(data['kla_1_lfc_trans'] < 1) | (data['dex_over_kla_1_lfc_trans'] > -.58)]
    up_in_kla = data[(data['kla_1_lfc_trans'] >= 1) & (data['dex_over_kla_1_lfc_trans'] > -.58)]
    
    
    
    
    for subgroup, suffix, dataset in (('Nearby Enhancers', '', data), 
                                      ('RefSeq Transcripts', '_trans', 
                                       data.groupby(by='nearest_refseq_transcript_id', as_index=False).mean())):
        
        supersets = (('All', dataset),
                 ('Not near transrepressed genes', not_trans), 
                 ('Up in KLA', up_in_kla),
                 ('Near transrepressed genes',transrepressed))
        
        
        labels = zip(*supersets)[0]
        xlabel = 'Subgroup'
    
        # KLA tag count
        vals = [d['kla_1_tag_count{0}'.format(suffix)]/d['length{0}'.format(suffix)] for _,d in supersets]
        title = 'GroSeq tag counts in {0}: DMSO 2h + KLA 1h by gene category'.format(subgroup)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='KLA 1h Tag Count per basepair', 
                         show_outliers=False, show_plot=False
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        yzer.show_plot()
        
        # KLA LFC
        vals = [d['kla_1_lfc{0}'.format(suffix)] for _,d in supersets]
        title = 'GroSeq Log Fold Change in {0} in DMSO 2h + KLA 1h by gene category'.format(subgroup)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='log2(KLA GRO-seq/DMSO GRO-seq)', 
                         show_outliers=False, show_plot=False
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        yzer.show_plot()
        
        # DMSO tag count
        vals = [d['dmso_1_tag_count{0}'.format(suffix)]/d['length{0}'.format(suffix)] for _,d in supersets]
        title = 'GroSeq tag counts in {0} in DMSO 2h by gene category'.format(subgroup)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='DMSO 2h Tag Count per basepair', 
                         show_outliers=False, show_plot=False
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        yzer.show_plot()
        
        # Dex LFC
        vals = [d['dex_1_lfc{0}'.format(suffix)] for _,d in supersets]
        title = 'GroSeq Log Fold Change in {0} in Dex 2h by gene category'.format(subgroup)
        ax = yzer.boxplot(vals, labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='log2(Dex GRO-seq/DMSO GRO-seq)', 
                         show_outliers=False, show_plot=False
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        yzer.show_plot()
        
        # p65 tag counts
        vals = [d['p65_kla_tag_count{0}'.format(suffix)] for _,d in supersets] \
                + [d['p65_kla_dex_tag_count{0}'.format(suffix)] for _,d in supersets]
        special_labels = ['{0}\np65 in KLA\nin {1}'.format(label, subgroup) for label,_ in supersets] \
                        + ['{0}\np65 in KLA+Dex\nin {1}'.format(label, subgroup) for label,_ in supersets]
        title = 'p65 Tag Counts in Peaks by gene category'
        ax = yzer.boxplot(vals, special_labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='Tag Counts in Peaks', 
                         show_outliers=False, show_plot=False, wide=True
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        yzer.show_plot()
        
        # GR tag counts
        vals = [d['gr_dex_tag_count{0}'.format(suffix)] for _,d in supersets] \
                + [d['gr_kla_dex_tag_count{0}'.format(suffix)] for _,d in supersets]
        special_labels = ['{0}\nGR in Dex\nin {1}'.format(label, subgroup) for label,_ in supersets] \
                        + ['{0}\nGR in KLA+Dex\nin {1}'.format(label, subgroup) for label,_ in supersets]
        title = 'GR Tag Counts in Peaks by gene category'
        ax = yzer.boxplot(vals, special_labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='Tag Counts in Peaks', 
                         show_outliers=False, show_plot=False, wide=True
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        yzer.show_plot()
        
        # PU.1 tag counts
        vals = [d['pu_1_tag_count{0}'.format(suffix)] for _,d in supersets] \
                + [d['pu_1_kla_tag_count{0}'.format(suffix)] for _,d in supersets]
        special_labels = ['{0}\nPU.1 in DMSO\nin {1}'.format(label, subgroup) for label,_ in supersets] \
                        + ['{0}\nPU.1 in KLA\nin {1}'.format(label, subgroup) for label,_ in supersets]
        title = 'PU.1 Tag Counts in Peaks by gene category'
        ax = yzer.boxplot(vals, special_labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='Tag Counts in Peaks', 
                         show_outliers=False, show_plot=False, wide=True
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        yzer.show_plot()