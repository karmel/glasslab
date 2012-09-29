'''
Created on Sep 26, 2012

@author: karmel
'''

from glasslab.dataanalysis.misc.gr_project_2012.width_buckets import get_data_with_bucket_score
from glasslab.dataanalysis.graphing.basepair_counter import BasepairCounter
from glasslab.dataanalysis.misc.gr_project_2012.elongation import set_up_sequencing_run_ids,\
    total_tags_per_run, get_rep_string

def draw_elongation_profile(data, grapher, dirpath, show_moving_average=True, show_count=False):
    run_ids = set_up_sequencing_run_ids()
    total_tags = total_tags_per_run()

    lfcs = (#('Special', 'group_{0}'),
              #('KLA','kla_{0}gene_body_lfc'),# ('KLA+Dex','kla_dex_{0}gene_body_lfc'),
              ('KLA+Dex over KLA','dex_over_kla_{0}gene_body_lfc'),
              )
    for desc,lfc in lfcs:
        for replicate_id in ('',1,2,3,4):
            rep_str = get_rep_string(replicate_id)
            lfc_str = lfc.format(rep_str)
            
            # Include all transcripts at once, but only do it once.
            if desc == 'Special': 
                datasets = [('All RefSeq', data),
                            ('Up > 2x in KLA, Down > 1.5x from that in Dex', 
                             data[(data['kla_{0}gene_body_lfc'.format(rep_str)] >= 1)
                                  & (data['dex_over_kla_{0}gene_body_lfc'.format(rep_str)] <= -1)]),]
            else:
                datasets = [('No change in {0}'.format(desc), data[data[lfc_str].abs() < 1]),
                            ('Up in {0}'.format(desc), data[data[lfc_str] >= 1]),
                            ('Down in {0}'.format(desc), data[data[lfc_str] <= -1]),]
            
            for label, dataset in datasets:
                slug_label = label.lower().replace(' ','_')
                group_by_cols = ['basepair','sequencing_run_id']
                data_grouped = dataset.groupby(group_by_cols, as_index=False).sum()
                
                groups = [data_grouped[data_grouped['sequencing_run_id'].isin(
                                                run_ids['dmso'][replicate_id or 0])],
                          data_grouped[data_grouped['sequencing_run_id'].isin(
                                                run_ids['kla'][replicate_id or 0])],
                          data_grouped[data_grouped['sequencing_run_id'].isin(
                                                run_ids['kla_dex'][replicate_id or 0])]]
                
                # Combine for sequencing runs now
                for i, group in enumerate(groups):
                    groups[i] = group.groupby(['basepair'], as_index=False).sum()
                    
                totals = zip(*total_tags.values())[replicate_id or 0]
                tag_scalars = grapher.get_tag_scalars(totals)
                ax = grapher.plot_tags_per_basepair(groups,
                            labels=['DMSO', 'KLA', 'KLA+Dex'],
                            title='Tag localization for RefSeq: {0}, {1}'.format(label,
                                    replicate_id and 'Group {0}'.format(replicate_id) or 'overall'),
                            tag_scalars=tag_scalars, show_moving_average=show_moving_average, 
                            show_count=show_count)
                grapher.save_plot(grapher.get_filename(dirpath, 
                            '{0}_refseq_by_run_type_{1}.png'.format(slug_label,
                                    replicate_id and 'group_{0}'.format(replicate_id) or 'all')))
                #grapher.show_plot()

if __name__ == '__main__':
    grapher = BasepairCounter()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/basepair_counts_gene_body'
    dirpath = grapher.get_path(dirpath)
    
    
    # Produce a dataframe with both basepair and overall gene body lfc for the transcript id
    if False:
        lfc_data = grapher.import_file(grapher.get_filename(dirpath,'feature_vectors.txt'))
        lfc_subset = lfc_data[['glass_transcript_id'] + filter(lambda x: '_gene_body_lfc' in x, lfc_data.columns)]
        bp_data = grapher.import_file(grapher.get_filename(dirpath,'refseq_by_basepair_and_run.txt'))
        data = bp_data.merge(lfc_subset, how='left', on='glass_transcript_id')
        data = data.fillna(0)
        data.to_csv(grapher.get_filename(dirpath,'refseq_by_basepair_with_gene_body_lfc.txt'),
                    sep='\t', header=True, index=False)
        
    data = grapher.import_file(grapher.get_filename(dirpath,'refseq_by_basepair_with_gene_body_lfc.txt'))
    draw_elongation_profile(data, grapher, dirpath)
    