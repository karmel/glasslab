'''
Created on Oct 17, 2012

@author: karmel


Do p65, PU.1, and GR all co-bind at these sites of interest?

Supershift to see if we can pick DNA locations, and pull out GR, p65, and PU.1 in KLA, KLA+Dex

Target: Has GR, p65, PU.1; loses p65, PU.1 with Dex
PosCtl: Has GR, p65, PU.1; no loss of p65, PU.1 with Dex
NegCtl: Has one of GR, p65, PU.1, but not the others

'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/'
    dirpath = yzer.get_path(dirpath)
    save_dirpath = yzer.get_and_create_path(dirpath,'subgroups_for_supershift','have_pu_1')
    
    transcripts = yzer.import_file(yzer.get_filename(dirpath, 'motifs', 'transcript_vectors.txt'))
    
    if True:
        data = transcripts[transcripts['refseq'] == 'f']
        data = data[data['has_infrastructure'] == 0]
        data = data[data['length'] < 6000]
        data = data[data['dex_1_lfc'] < 1]
        data = data[data['kla_1_lfc'] >= 1]
        
        data = data.fillna(0)
        
        # Redo our links to be more useful for display
        data['ucsc_link_nod'] = data['ucsc_link_nod'].map(lambda x: '<a href={0} target="_blank">UCSC</a>'.format(
                                                x.replace('nod_balbc','gr_project_2012')))
        
        # First get sets for Negative controls
        tfs = ['p65','pu_1','gr','gr_fa']
        for tf in tfs:
            other_tfs = tfs[:]
            other_tfs.remove(tf)
        
            # Get only TF of interest
            only = data[data['{0}_kla_dex_tag_count'.format(tf)] > 0]
            for other in other_tfs:
                only = only[only['{0}_kla_dex_tag_count'.format(other)] == 0]
        
            only.to_csv(yzer.get_filename(save_dirpath, 'enhancer_like_{0}_only.txt'.format(tf)), 
                                      sep='\t', header=True, index=False)
            
        
        # Next look only at sites where we have all three binding, but not GR directly.
        data = data[data['gr_kla_dex_tag_count'] > 0]
        #data = data[data['pu_1_kla_dex_tag_count'] + data['pu_1_kla_tag_count'] > 0]
        data = data[data['p65_kla_dex_tag_count'] + data['p65_kla_tag_count'] > 0]
        data = data[data['pu_1_kla_dex_tag_count'] + data['pu_1_kla_tag_count'] > 0]
        data = data[data['gr_fa_kla_dex_tag_count'] == 0]
        
        # Split up by change in p65
        peak_type = 'p65'
        peak_type_2 = 'pu_1'
        ratio = 3
        ratio_2 = 1.2
        kla_gt = (data['{0}_kla_tag_count'.format(peak_type)] > ratio*data['{0}_kla_dex_tag_count'.format(peak_type)])
        kla_gt_both = kla_gt & (data['{0}_kla_tag_count'.format(peak_type_2)] > ratio_2*data['{0}_kla_dex_tag_count'.format(peak_type_2)])
        
        nc = (data['{0}_kla_tag_count'.format(peak_type)] + data['{0}_kla_dex_tag_count'.format(peak_type)] > 0) \
            & (data['{0}_kla_dex_tag_count'.format(peak_type)] < ratio*data['{0}_kla_tag_count'.format(peak_type)]) \
            & (data['{0}_kla_tag_count'.format(peak_type)] < ratio*data['{0}_kla_dex_tag_count'.format(peak_type)])
        
        nc_both = nc & (data['{0}_kla_tag_count'.format(peak_type_2)] + data['{0}_kla_dex_tag_count'.format(peak_type_2)] > 0) \
            & (data['{0}_kla_dex_tag_count'.format(peak_type_2)] < ratio_2*data['{0}_kla_tag_count'.format(peak_type_2)]) \
            & (data['{0}_kla_tag_count'.format(peak_type_2)] < ratio_2*data['{0}_kla_dex_tag_count'.format(peak_type_2)])
        
        data['{0}_kla_to_kla_dex_ratio'.format(peak_type)] = data['{0}_kla_tag_count'.format(peak_type)]\
                                                                /data['{0}_kla_dex_tag_count'.format(peak_type)]
        
        data[kla_gt].to_csv(yzer.get_filename(save_dirpath, 'enhancer_like_lose_{0}_{1}x_change_dsg_only.txt'.format(
                                                                                        peak_type, ratio)), 
                                      sep='\t', header=True, index=False)
        data[kla_gt_both].to_csv(yzer.get_filename(save_dirpath, 'enhancer_like_lose_{0}_{2}_{1}x_change_dsg_only.txt'.format(
                                                                                        peak_type, ratio, peak_type_2)), 
                                      sep='\t', header=True, index=False)
        data[nc].to_csv(yzer.get_filename(save_dirpath, 'enhancer_like_{0}_no_change_dsg_only.txt'.format(peak_type)), 
                                      sep='\t', header=True, index=False)
        data[nc_both].to_csv(yzer.get_filename(save_dirpath, 'enhancer_like_{0}_{1}_no_change_dsg_only.txt'.format(
                                                                                        peak_type, peak_type_2)), 
                                      sep='\t', header=True, index=False)
        