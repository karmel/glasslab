'''
Created on Mar 23, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/NOD_BALBc/ThioMacs/Analysis_2013_02/'
    dirpath = yzer.get_path(dirpath)
    img_dirpath = yzer.get_and_create_path(dirpath, 'refseq_expression')
    data = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
    
    data = data.fillna(0)
    data = yzer.normalize(data, 'nod_notx_1h_tag_count', 1.095436)
    data = yzer.normalize(data, 'nod_kla_1h_tag_count', 0.652898)
    #data = yzer.normalize(data, 'nonplated_diabetic_nod_notx_tag_count', 0.885427)
    #data = yzer.normalize(data, 'nonplated_diabetic_balb_notx_tag_count', 0.645579)
    
    data['balb_notx_1h_reads_per_base'] = data['balb_notx_1h_tag_count']/data['length']
    data['balb_kla_1h_reads_per_base'] = data['balb_kla_1h_tag_count']/data['length']
    
    refseq = yzer.get_refseq(data)
    
    # Remove low tag counts
    refseq = refseq[refseq['transcript_score'] >= 4]
    

    if False:
        # Non-diabetic balbc vs. nod
        ax = yzer.scatterplot(refseq, 'balb_notx_1h_tag_count', 'nod_notx_1h_tag_count_norm',
                            log=True, color='blue', master_dataset=refseq,
                            title='Non-Diabetic BALBc vs. NOD Notx 1h Refseq Transcripts',
                            show_2x_range=True, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)
        
        # Non-diabetic balbc vs. nod
        ax = yzer.scatterplot(refseq, 'balb_kla_1h_tag_count', 'nod_kla_1h_tag_count_norm',
                            log=True, color='blue', master_dataset=refseq,
                            title='Non-Diabetic BALBc vs. NOD KLA 1h Refseq Transcripts',
                            show_2x_range=True, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)
        
    
            
    if True:
        gene_groups = [['Clec4e','Tlr2',],
                       ['Cxcl1','Cxcl2','Il6','Ptgs2','Tnfsf9','Vegfa','Tnf',
                 'Siglec1','Mmp9',
                 'Il10','Il1b','Cxcl10','Tlr4','Il12b',]]
        '''
                 'Itgb2','Itgam','Rhoa','Limk1',
                 'Cdc42','Rac1','Rac2','Vav1','Vav2','Vav3',
                 'Arpc1a','Arpc1b','Arpc2','Arpc3','Arpc5','Actr2','Actr3',
                 'Cfl1','Cfl2','Dnm2','Arf1','Arf6',
                 'Pdlim5', 'Pls3', 'Lima1', 'Vcl', 'Ctnnb1', 'Coro1c', 'Cdc42bpg', 'Flna', 'Flnc', 
                 'Cnn2', 'Coro1a', 'Myl6', 'Cfl1', 'Csrp1', 'Srf', 'Cnn3', 'Arhgef17', 'Myh9', 'Myo18b', 'Mybpc3']
        '''
        for i, genes in enumerate(gene_groups):
            indices = [refseq[refseq['gene_names'] == gene].index[0] for gene in genes]
            
            
            sorted_by_count = refseq.fillna(0).sort_index(axis=0, by='balb_kla_1h_reads_per_base').index.copy()
            sort_indexes = list(enumerate(sorted_by_count))
            sort_indexes.sort(key=lambda x: x[1])
            refseq['rank'] = zip(*sort_indexes)[0]
             
            yzer.bargraph_for_transcripts(refseq, indices, ['balb_nod_notx_1h_fc'],
                                                bar_names=genes,
                                                title='Fold Change in NOD vs. BALBc notx 1h GRO-seq',
                                                ylabel='Fold Change in NOD vs. BALBc',
                                                rank_label='Rank of read per base pair value in BALBc notx 1h, ascending',
                                                show_plot=False)
            yzer.save_plot(yzer.get_filename(img_dirpath, 'balbc_nod_notx_' + str(i) + '_fold_change_bargraph.png'))
            yzer.show_plot()
                