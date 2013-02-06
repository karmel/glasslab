'''
Created on Jul 2, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.basepair_counter import BasepairCounter
import os
from glasslab.dataanalysis.misc.gr_project_2012.v1.elongation import set_up_sequencing_run_ids,\
    total_tags_per_run, get_rep_string


def draw_elongation_profile(data, grapher, dirpath, show_moving_average=True, show_count=False):
    run_ids = set_up_sequencing_run_ids()
    total_tags = total_tags_per_run()

    states = (('Special', 'group_{0}'),
              ('KLA','kla_{0}lfc'), ('KLA+Dex','kla_dex_{0}lfc'),
              ('KLA+Dex over KLA','dex_over_kla_{0}lfc'),
              )
    for desc,state in states:
        for replicate_id in ('',1,2,3,4):
            rep_str = get_rep_string(replicate_id)
            state_str = state.format(rep_str)
            
            # Include all transcripts at once, but only do it once.
            if desc == 'Special': 
                datasets = [('All RefSeq', data),
                            ('Up > 2x in KLA, Down > 1.5x from that in Dex', 
                             data[(data['kla_{0}lfc'.format(rep_str)] >= 1)
                                  & (data['dex_over_kla_{0}lfc'.format(rep_str)] <= -.58)]),]
            else:
                datasets = [('No change in {0}'.format(desc), data[data[state_str].abs() < 1]),
                            ('Up in {0}'.format(desc), data[data[state_str] >= 1]),
                            ('Down in {0}'.format(desc), data[data[state_str] <= -1]),]
            
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
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/buckets_per_row_low_basal'
    dirpath = grapher.get_path(dirpath)
    
    if True: 
        data = grapher.import_file(grapher.get_filename(dirpath, 'transcripts_per_run_per_bucket.txt'))
        # For the sake of graphing, imitate basepair
        data['basepair'] = (data['bucket'] - 1)*50

        # Merge on transcripts to get DMSO tags.
        transcripts = grapher.import_file(grapher.get_filename(dirpath, 'transcript_vectors.txt'))
        transcripts = transcripts[transcripts['has_refseq'] == 1]
        transcripts['glass_transcript_id'] = transcripts['id']
        totals = total_tags_per_run()
        for replicate_id in ('',1,2,3,4):
            rep_str = get_rep_string(replicate_id)
            transcripts['dmso_{0}rpkm'.format(rep_str)] = transcripts['dmso_{0}tag_count'.format(rep_str)]*(10**3*10**6)\
                        /transcripts['length']/totals['dmso'][replicate_id or 0]    
                    
        transcripts = transcripts.filter(regex=r'glass_transcript_id|has_refseq|.*lfc|.*rpkm')
        
        data = data.merge(transcripts, how='left', on='glass_transcript_id')
        data = data.fillna(0)
        data = data[data['has_refseq'] == 1]
        
        # Create filtered groups.
        datasets = [('all_refseq', data),
                    ('low_basal', data[data['dmso_rpkm'] < 1]),
                    ('low_basal_1', data[data['dmso_1_rpkm'] < 1]),
                    ('low_basal_2', data[data['dmso_2_rpkm'] < 1]),
                    ('low_basal_3', data[data['dmso_3_rpkm'] < 1]),
                    ('low_basal_4', data[data['dmso_4_rpkm'] < 1]),
                    ]
    
    for name, dataset in datasets:
            
     
            curr_dirpath = grapher.get_filename(dirpath, name)
            if not os.path.exists(curr_dirpath): os.mkdir(curr_dirpath)
            
            draw_elongation_profile(dataset, grapher, curr_dirpath, 
                                    show_moving_average=False, show_count=True)

            
            