'''
Created on Feb 14, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/Oshea_enhancers/ctcf_across_celltypes'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'figures')
    
    dp = yzer.import_file(yzer.get_filename(dirpath, 'dp_with_thiomac_ctcf.txt')).fillna(0)
    thio = yzer.import_file(yzer.get_filename(dirpath, 'thiomac_with_dp_ctcf.txt')).fillna(0)
    
    # Get venn-diagram sets
    only_dp = dp[dp['thiomac_ctcf_tag_count'] == 0]
    only_thio = thio[thio['dp_ctcf_tag_count'] == 0]
    shared = dp[dp['thiomac_ctcf_tag_count'] != 0]
    shared_check = thio[thio['dp_ctcf_tag_count'] != 0]
    print len(only_dp), len(only_thio), len(shared), len(shared_check)
    
    data = shared.append(only_dp, ignore_index=True)
    data = data.append(only_thio, ignore_index=True)

    data['dp_nonzero'] = nonzero(data['dp_ctcf_tag_count'])
    data['thio_nonzero'] = nonzero(data['thiomac_ctcf_tag_count'])
    ax = yzer.scatterplot(data, 'dp_nonzero', 'thio_nonzero',
                            xlabel='DP Thymocyte CTCF Tag Count', ylabel='ThioMac CTCF Tag Count',
                            log=True, color='blue', 
                            title='Tags in CTCF Peaks in DP Thymocytes versus ThioMacs',
                            show_2x_range=False, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)