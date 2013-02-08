'''
Created on Oct 26, 2012

@author: karmel
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from matplotlib import pyplot
from glasslab.utils.functions import nonzero
from glasslab.dataanalysis.misc.gr_project_2012.v1.elongation import total_tags_per_run

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/enhancer_classification'
    dirpath = yzer.get_path(dirpath)
    
    consistent = False
    img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_by_expression', consistent and 'consistent' or 'rep1')
    
    data = yzer.import_file(yzer.get_filename(dirpath, 'enhancers_with_nearest_gene.txt'))
    data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
    
    draw_pies = True
    min_tags = 30
    ratio = 1.5
    # Make sure we have dimethyl
    data = data[data.filter(like='h3k4me2').max(axis=1) > min_tags]
    data = data[data['minimal_distance'] >= 1000]
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    transcripts['nearest_refseq_transcript_id'] = transcripts['id']
    data = data.merge(transcripts, how='left', on='nearest_refseq_transcript_id', suffixes=['','_trans'])
    
    data = data.fillna(0)
    
    total_tags = total_tags_per_run()
        
    data['h4k8ac_kla_dex_ratio'] = nonzero(data['h4k8ac_kla_dex_tag_count'])/nonzero(data['h4k8ac_kla_tag_count'])
    data['dmso_1_rpkm'] = data['dmso_1_tag_count_trans']*(10**3*10**6)/data['length_trans']/total_tags['dmso'][1]
        
    if False:
        # Low basal expression only
        data = data[data['dmso_1_rpkm'] > 2]
        img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_by_expression_high_basal', consistent and 'consistent' or 'rep1')
    if False:
        # Lose acetyl only
        data = data[data['h4k8ac_kla_dex_ratio'] < .75]
        img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_by_expression_lose_ac', consistent and 'consistent' or 'rep1')
        
    
    
    
    
    
    for subgroup, suffix, dataset in (('Nearby Enhancers', '', data), 
                                      ('RefSeq Transcripts', '_trans', 
                                       data.groupby(by='nearest_refseq_transcript_id', as_index=False).mean())):
        if consistent:
            transrepressed = dataset[(dataset['kla_1_lfc_trans'] >= 1) & (dataset['dex_over_kla_1_lfc_trans'] <= -.58)
                                     & (dataset['kla_3_lfc_trans'] >= 1) & (dataset['dex_over_kla_3_lfc_trans'] <= -.58)
                                     & (dataset['kla_4_lfc_trans'] >= 1) & (dataset['dex_over_kla_4_lfc_trans'] <= -.58)]
            up_in_kla = dataset[((dataset['kla_1_lfc_trans'] >= 1) & (dataset['kla_3_lfc_trans'] >= 1) 
                                & (dataset['kla_4_lfc_trans'] >= 1)) & 
                                ((dataset['dex_over_kla_1_lfc_trans'] > -.58) | (dataset['dex_over_kla_3_lfc_trans'] > -.58)
                                 | (dataset['dex_over_kla_4_lfc_trans'] > -.58))]
        else:
            transrepressed = dataset[(dataset['kla_1_lfc_trans'] >= 1) & (dataset['dex_over_kla_1_lfc_trans'] <= -.58)]
            up_in_kla = dataset[(dataset['kla_1_lfc_trans'] >= 1) & (dataset['dex_over_kla_1_lfc_trans'] > -.58)]
        
        not_trans = dataset[(~dataset.index.isin(up_in_kla.index)) & (~dataset.index.isin(transrepressed.index))]
        
        supersets = (('All', dataset),
                 ('Not {0}Up in KLA Genes'.format(not suffix and 'Near ' or ''), not_trans), 
                 ('{0}Up in KLA Genes'.format(not suffix and 'Near ' or ''), up_in_kla),
                 ('{0}Transrepressed Genes'.format(not suffix and 'Near ' or ''),transrepressed))
        
        print len(up_in_kla), len(transrepressed), len(not_trans)
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
        
        
        # p65 tag counts
        vals = [d['p65_kla_tag_count{0}'.format(suffix)] for _,d in supersets] \
                + [d['p65_kla_dex_tag_count{0}'.format(suffix)] for _,d in supersets]
        special_labels = ['{0}\np65 in KLA\nin {1}'.format(label, subgroup) for label,_ in supersets] \
                        + ['{0}\np65 in KLA+Dex\nin {1}'.format(label, subgroup) for label,_ in supersets]
        title = 'p65 Tag Counts in Peaks in {0} by gene category'.format(subgroup)
        ax = yzer.boxplot(vals, special_labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='Tag Counts in Peaks', 
                         show_outliers=False, show_plot=False, wide=True
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        
        
        # GR tag counts
        vals = [d['gr_dex_tag_count{0}'.format(suffix)] for _,d in supersets] \
                + [d['gr_kla_dex_tag_count{0}'.format(suffix)] for _,d in supersets]
        special_labels = ['{0}\nGR in Dex\nin {1}'.format(label, subgroup) for label,_ in supersets] \
                        + ['{0}\nGR in KLA+Dex\nin {1}'.format(label, subgroup) for label,_ in supersets]
        title = 'GR Tag Counts in Peaks in {0} by gene category'.format(subgroup)
        ax = yzer.boxplot(vals, special_labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='Tag Counts in Peaks', 
                         show_outliers=False, show_plot=False, wide=True
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        
        
        # PU.1 tag counts
        vals = [d['pu_1_tag_count{0}'.format(suffix)] for _,d in supersets] \
                + [d['pu_1_kla_tag_count{0}'.format(suffix)] for _,d in supersets]
        special_labels = ['{0}\nPU.1 in DMSO\nin {1}'.format(label, subgroup) for label,_ in supersets] \
                        + ['{0}\nPU.1 in KLA\nin {1}'.format(label, subgroup) for label,_ in supersets]
        title = 'PU.1 Tag Counts in Peaks in {0} by gene category'.format(subgroup)
        ax = yzer.boxplot(vals, special_labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='Tag Counts in Peaks', 
                         show_outliers=False, show_plot=False, wide=True
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        
        
        # H3K4me2 tag counts
        vals = [d['h3k4me2_tag_count{0}'.format(suffix)] for _,d in supersets] \
                + [d['h3k4me2_kla_tag_count{0}'.format(suffix)] for _,d in supersets]
        special_labels = ['{0}\nH3K4me2 in DMSO\nin {1}'.format(label, subgroup) for label,_ in supersets] \
                        + ['{0}\nH3K4me2 in KLA\nin {1}'.format(label, subgroup) for label,_ in supersets]
        title = 'H3K4me2 Tag Counts in Peaks in {0} by gene category'.format(subgroup)
        ax = yzer.boxplot(vals, special_labels, 
                         title=title, xlabel=xlabel, 
                         ylabel='Tag Counts in Peaks', 
                         show_outliers=False, show_plot=False, wide=True
                         )
        pyplot.setp(ax.get_xticklabels(), fontsize=10)
        yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
        
        
    # H4K8ac tag counts
    vals = [d['h4k8ac_tag_count'.format(suffix)] for _,d in supersets] \
            + [d['h4k8ac_kla_tag_count'.format(suffix)] for _,d in supersets]
    special_labels = ['{0}\nH4K8ac in DMSO\nin {1}'.format(label, subgroup) for label,_ in supersets] \
                    + ['{0}\nH4K8ac in KLA\nin {1}'.format(label, subgroup) for label,_ in supersets]
    title = 'H4K8ac Tag Counts in Peaks in {0} by gene category'.format(subgroup)
    ax = yzer.boxplot(vals, special_labels, 
                     title=title, xlabel=xlabel, 
                     ylabel='Tag Counts in Peaks', 
                     show_outliers=False, show_plot=False, wide=True
                     )
    pyplot.setp(ax.get_xticklabels(), fontsize=10)
    yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
    
    
    # H4K8ac tag changes
    vals = [nonzero(d['h4k8ac_kla_tag_count'])/nonzero(d['h4k8ac_tag_count']) for _,d in supersets]
    special_labels = ['{0}\nH4K8ac Ratio\nin {1}'.format(label, subgroup) for label,_ in supersets] 
    title = 'H4K8ac Tag Count Ratios in Peaks in {0} by gene category'.format(subgroup)
    ax = yzer.boxplot(vals, special_labels, 
                     title=title, xlabel=xlabel, 
                     ylabel='KLA Tags/DMSO Tags', 
                     show_outliers=False, show_plot=False, wide=False
                     )
    pyplot.setp(ax.get_xticklabels(), fontsize=10)
    yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
    
    # H4K8ac tag counts
    vals = [d['h4k8ac_tag_count'.format(suffix)] for _,d in supersets] \
            + [d['h4k8ac_kla_dex_tag_count'.format(suffix)] for _,d in supersets]
    special_labels = ['{0}\nH4K8ac in DMSO\nin {1}'.format(label, subgroup) for label,_ in supersets] \
                    + ['{0}\nH4K8ac in Dex+KLA\nin {1}'.format(label, subgroup) for label,_ in supersets]
    title = 'Dex+KLA H4K8ac Tag Counts in Peaks in {0} by gene category'.format(subgroup)
    ax = yzer.boxplot(vals, special_labels, 
                     title=title, xlabel=xlabel, 
                     ylabel='Tag Counts in Peaks', 
                     show_outliers=False, show_plot=False, wide=True
                     )
    pyplot.setp(ax.get_xticklabels(), fontsize=10)
    yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
    
    
    # H4K8ac tag changes
    vals = [nonzero(d['h4k8ac_kla_dex_tag_count'])/nonzero(d['h4k8ac_tag_count']) for _,d in supersets]
    special_labels = ['{0}\nH4K8ac Ratio\nin {1}'.format(label, subgroup) for label,_ in supersets] 
    title = 'Dex+KLA H4K8ac Tag Count Ratios in Peaks in {0} by gene category'.format(subgroup)
    ax = yzer.boxplot(vals, special_labels, 
                     title=title, xlabel=xlabel, 
                     ylabel='Dex+KLA Tags/DMSO Tags', 
                     show_outliers=False, show_plot=False, wide=False
                     )
    pyplot.setp(ax.get_xticklabels(), fontsize=10)
    yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
    
    # H4K8ac tag changes
    vals = [nonzero(d['h4k8ac_kla_dex_tag_count'])/nonzero(d['h4k8ac_kla_tag_count']) for _,d in supersets]
    special_labels = ['{0}\nH4K8ac Ratio\nin {1}'.format(label, subgroup) for label,_ in supersets] 
    title = 'Dex+KLA to KLA H4K8ac Tag Count Ratios in Peaks in {0} by gene category'.format(subgroup)
    ax = yzer.boxplot(vals, special_labels, 
                     title=title, xlabel=xlabel, 
                     ylabel='Dex+KLA Tags/KLA Tags', 
                     show_outliers=False, show_plot=False, wide=False
                     )
    pyplot.setp(ax.get_xticklabels(), fontsize=10)
    yzer.save_plot(yzer.get_filename(img_dirpath, title + '.png'))
    