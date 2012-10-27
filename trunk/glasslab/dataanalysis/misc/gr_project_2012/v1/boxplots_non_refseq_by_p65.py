'''
Created on Oct 1, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero
from glasslab.dataanalysis.misc.gr_project_2012.enhancer_subsets_for_supershift import ucsc_link_cleanup
import numpy

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    peak_type = 'p65'
        
    img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_non_refseq_by_{0}'.format(peak_type))
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'transcript_vectors_with_nearby_peaks.txt'))
    
    
    if True:
        pu_1 = False
        for ratio in (1.5, 2, 3):
            data = transcripts[transcripts['refseq'] == 'f']
            data = data[data['has_infrastructure'] == 0]
            data = data[data['length'] < 6000]
            data = data[data['dex_1_lfc'] < 1]
            data = data[data['kla_1_lfc'] >= 1]
            data = data[data['gr_kla_dex_tag_count'] > 0]
            data = data[data['gr_fa_kla_dex_tag_count'] == 0]
            print len(data)
            if pu_1: data = data[data['pu_1_kla_tag_count']  + data['pu_1_kla_tag_count'] > 0]
            
            
            data = data.fillna(0)
            data = ucsc_link_cleanup(data)
            
            
            colname = 'dex_over_kla_1_lfc'
            
            if pu_1:
                colname = 'pu_1_ratio'
                data[colname] = numpy.log2(nonzero(data['pu_1_kla_dex_tag_count'])/nonzero(data['pu_1_kla_tag_count']))
            
            none = (data['{0}_kla_nearby_tag_count'.format(peak_type)] + data['{0}_kla_dex_nearby_tag_count'.format(peak_type)] == 0)
            kla_gt = (data['{0}_kla_tag_count'.format(peak_type)] > ratio*data['{0}_kla_dex_nearby_tag_count'.format(peak_type)])
            kla_dex_gt = (data['{0}_kla_dex_tag_count'.format(peak_type)] > ratio*data['{0}_kla_nearby_tag_count'.format(peak_type)])
            nc = (data['{0}_kla_nearby_tag_count'.format(peak_type)] + data['{0}_kla_dex_nearby_tag_count'.format(peak_type)] > 0) \
                & (data['{0}_kla_dex_tag_count'.format(peak_type)] < ratio*data['{0}_kla_nearby_tag_count'.format(peak_type)]) \
                & (data['{0}_kla_tag_count'.format(peak_type)] < ratio*data['{0}_kla_dex_nearby_tag_count'.format(peak_type)])
            
            data['{0}_kla_to_kla_dex_ratio'.format(peak_type)] = data['{0}_kla_tag_count'.format(peak_type)]\
                                                                    /data['{0}_kla_dex_nearby_tag_count'.format(peak_type)]
            
            print sum(kla_gt)
            print sum(kla_dex_gt)
            data[kla_gt].to_csv(yzer.get_filename(img_dirpath, 'enhancer_like_lose_{0}_{1}x_change_dsg_only.txt'.format(
                                                                                            peak_type, ratio)), 
                                          sep='\t', header=True, index=False)
            data[kla_dex_gt].to_csv(yzer.get_filename(img_dirpath, 'enhancer_like_gain_{0}_{1}x_change_dsg_only.txt'.format(
                                                                                            peak_type, ratio)), 
                                          sep='\t', header=True, index=False)
            
            title = 'LFC in KLA + Dex over KLA by change in {0}:\nEnhancer-Like, Has GR in KLA+Dex (DSG only)'.format(peak_type)
            if pu_1:
                title = 'PU.1 in KLA + Dex over KLA by change in {0}:\nEnhancer-Like, has GR in KLA+Dex, has PU.1'.format(peak_type)
            
            names = [s.format(peak_type) for s in ['No {0}','Loses {0} in KLA+Dex','No change in {0}', 'Gains {0} in KLA+Dex']]
            ax = yzer.boxplot([data[none][colname], data[kla_gt][colname], data[nc][colname], data[kla_dex_gt][colname]], 
                         names,
                         title=title, 
                         xlabel='{0} Status'.format(peak_type), 
                         ylabel=(pu_1 and 'log2(KLA+Dex PU.1/KLA PU.1)')\
                            or 'log2(KLA+Dex GRO-seq/KLA GRO-seq)', 
                         show_outliers=False, show_plot=False)
            if pu_1:
                yzer.save_plot(yzer.get_filename(img_dirpath, 'dex_over_kla_pu_1_both_by_{0}_sums_{1}x_change.png'.format(peak_type, ratio)))
            else:
                yzer.save_plot(yzer.get_filename(img_dirpath, 'dex_over_kla_1_lfc_by_{0}_sums_{1}x_change_dsg_only.png'.format(peak_type, ratio)))
            
            yzer.show_plot()