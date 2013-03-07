'''
Created on Mar 4, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/NOD_BALBc/ThioMacs/Analysis_2013_02/'
    dirpath_bmdc = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/NOD_BALBc/BMDCs/Analysis_2013_03/'
    dirpath = yzer.get_path(dirpath)
    dirpath_bmdc = yzer.get_path(dirpath_bmdc)
    img_dirpath = yzer.get_and_create_path(dirpath, 'bmdc_vs_thiomac')
    thio = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    bmdc = yzer.import_file(yzer.get_filename(dirpath_bmdc, 'transcript_vectors.txt'))
    
    sets = []
    
    for data in (thio, bmdc):
        data = data.fillna(0)
    
        refseq = yzer.get_refseq(data)
        
        # Remove low tag counts
        #refseq = refseq[refseq['transcript_score'] >= 4]
        
        sets.append(refseq)
        
    
            
    if True:
        genes = ['Coro1a','Vcl','Tlr2','Clec4e','Cxcl2']
        
        vals = []
        labels = []
        for gene in genes: 
            vals.append(sets[0][sets[0]['gene_names'] == gene]['balb_nod_notx_1h_fc'].values[0])

            vals.append(sets[1][sets[1]['gene_names'] == '{' + gene + '}']['balb_nod_notx_0h_fc'].values[0])
            labels.append('ThioMac\n' + gene)
            labels.append('BMDC\n' + gene)
        
        # Convert to FC from log
        vals = [2**val for val in vals]
        yzer.bargraph_for_transcripts(vals, None, None,
                        bar_names=labels,
                        title='Fold Change in NOD vs. BALBc NOTX GRO-seq: ThioMac and BMDC',
                        ylabel='Fold Change in NOD vs. BALBc',
                        show_plot=True, save_dir=img_dirpath)
        