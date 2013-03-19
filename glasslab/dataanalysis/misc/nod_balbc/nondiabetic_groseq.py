'''
Created on Mar 23, 2012

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
from glasslab.utils.functions import nonzero

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
    
    data['balb_notx_1h_tag_count'] = nonzero(data['balb_notx_1h_tag_count'])
    data['nod_notx_1h_tag_count_norm'] = nonzero(data['nod_notx_1h_tag_count_norm'])
    data['balb_kla_1h_tag_count'] = nonzero(data['balb_kla_1h_tag_count'])
    data['nod_kla_1h_tag_count_norm'] = nonzero(data['nod_kla_1h_tag_count_norm'])
    
    data = data[data['transcript_score'] >= 4]
    data = data[data[['balb_notx_1h_tag_count','nod_notx_1h_tag_count_norm',
                      'balb_kla_1h_tag_count','nod_kla_1h_tag_count_norm']].max(axis=1) >= 10]
    
    refseq = yzer.get_refseq(data)
    
    # Remove low tag counts
    refseq = refseq[refseq['transcript_score'] >= 4]
    

    if False:
        # Non-diabetic balbc vs. nod
        ax = yzer.scatterplot(refseq, 'balb_notx_1h_tag_count', 'nod_notx_1h_tag_count_norm',
                            log=True, color='blue', master_dataset=refseq,
                            xlabel='BALBc notx 1h tag count',ylabel='NOD notx 1h tag count',
                            title='Non-Diabetic BALBc vs. NOD Notx 1h Refseq Transcripts',
                            show_2x_range=True, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)
        
        # Non-diabetic balbc vs. nod
        ax = yzer.scatterplot(refseq, 'balb_kla_1h_tag_count', 'nod_kla_1h_tag_count_norm',
                            log=True, color='blue', master_dataset=refseq,
                            xlabel='BALBc KLA 1h tag count',ylabel='NOD KLA 1h tag count',
                            title='Non-Diabetic BALBc vs. NOD KLA 1h Refseq Transcripts',
                            show_2x_range=True, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)
    if False:
        # Non-diabetic balbc vs. nod
        ax = yzer.scatterplot(data, 'balb_notx_1h_tag_count', 'nod_notx_1h_tag_count_norm',
                            log=True, color='blue', master_dataset=data,
                            xlabel='BALBc notx 1h tag count',ylabel='NOD notx 1h tag count',
                            title='Non-Diabetic BALBc vs. NOD Notx 1h All Transcripts',
                            show_2x_range=True, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)
        
        # Non-diabetic balbc vs. nod
        ax = yzer.scatterplot(data, 'balb_kla_1h_tag_count', 'nod_kla_1h_tag_count_norm',
                            log=True, color='blue', master_dataset=data,
                            xlabel='BALBc KLA 1h tag count',ylabel='NOD KLA 1h tag count',
                            title='Non-Diabetic BALBc vs. NOD KLA 1h All Transcripts',
                            show_2x_range=True, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)
        
    if False:
        # Enhancers
        enh = data[data['distal'] == 't']
        enh = enh[enh['h3k4me2_tag_count'] > 10]
        # Non-diabetic balbc vs. nod
        ax = yzer.scatterplot(enh, 'balb_notx_1h_tag_count', 'nod_notx_1h_tag_count_norm',
                            log=True, color='blue', master_dataset=enh,
                            xlabel='BALBc notx 1h tag count',ylabel='NOD notx 1h tag count',
                            title='Non-Diabetic BALBc vs. NOD Notx 1h eRNA',
                            show_2x_range=True, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)
        
        # Non-diabetic balbc vs. nod
        ax = yzer.scatterplot(enh, 'balb_kla_1h_tag_count', 'nod_kla_1h_tag_count_norm',
                            log=True, color='blue', master_dataset=enh,
                            xlabel='BALBc KLA 1h tag count',ylabel='NOD KLA 1h tag count',
                            title='Non-Diabetic BALBc vs. NOD KLA 1h eRNA',
                            show_2x_range=True, show_legend=False,
                            show_count=True, show_correlation=True, 
                            save_dir=img_dirpath, show_plot=True)
        
    
            
    if True:
        gene_groups = {
                    #'tlr2_targets': ['Tlr2', 'Il10','Cxcl10','Il12b','Ccl5','Tlr4'],
                    'clec4e_targets_notx': ['Cxcl1', 'Cxcl2','Ptgs2','Tnf','Il1b'],
                    'clec4e_targets_kla': ['Cxcl1', 'Cxcl2','Ptgs2','Tnf','Il1b','Il6','Il12b',''],
                       }         
        '''
        'srf_targets': ['Srf','Cnn2','Lima1','Coro1a','Vcl','Acta2','Actb','Dhcr24','Actg2','Actc1','Lcp1','Jup','Tpm4','Tnni2','Zyx','Tubb3','Pfn1','Gas7','Arpc4','Pstpip1','Bsn','Flna','Actn1'],
                       
                'clec4e_tlr2': ['Clec4e','Tlr2',],
                       'inflammatory_genes': ['Cxcl1','Cxcl2','Il6','Ptgs2','Tnfsf9','Vegfa','Tnf',
                 'Siglec1','Mmp9',
                 'Il10','Il1b','Cxcl10','Tlr4','Il12b',],
        '''
        for i, genes in gene_groups.items():
            for gene in genes[:]:
                if not gene in refseq['gene_names'].values: 
                    print gene
                    genes.remove(gene)
            indices = [refseq[refseq['gene_names'] == gene].index[0] for gene in genes]
            
            for txt in ('notx','kla'):
                sorted_by_count = refseq.fillna(0).sort_index(axis=0, by='balb_{0}_1h_reads_per_base'.format(txt)).index.copy()
                sort_indexes = list(enumerate(sorted_by_count))
                sort_indexes.sort(key=lambda x: x[1])
                refseq['rank'] = zip(*sort_indexes)[0]
                 
                yzer.bargraph_for_transcripts(refseq, indices, ['balb_nod_{0}_1h_fc'.format(txt)],
                                bar_names=genes,
                                title='Fold Change in NOD vs. BALBc {0} 1h GRO-seq'.format(txt.upper()),
                                ylabel='Fold Change in NOD vs. BALBc',
                                rank_label='Rank of read per base pair value in BALBc {0} 1h, ascending'.format(txt.upper()),
                                show_plot=False)
                yzer.save_plot(yzer.get_filename(img_dirpath, 'balbc_nod_{0}_{1}_fold_change_bargraph.png'.format(txt,i)))
                #yzer.show_plot()
                