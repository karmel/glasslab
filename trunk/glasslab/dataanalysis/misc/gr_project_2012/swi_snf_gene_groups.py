'''
Created on Nov 26, 2012

@author: karmel

Do the gene groups outlined in Ramirez-Carrozzi 2006 and 2009 correlate
with expression changes in Dex+KLA?
'''
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher


if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/GR_Analysis/cpg_island_promoters'
    dirpath = yzer.get_path(dirpath)
    
    for rep in (1,3,4):
        img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_by_expression', 'Ramirez-Carrozzi','rep{0}'.format(rep))
        
        data = yzer.import_file(yzer.get_filename(dirpath, 'transcript_vectors.txt'))
        data['ucsc_link_nod'] = data['ucsc_link_nod'].apply(lambda s: s.replace('nod_balbc','gr_project_2012'))
        data = data.fillna(0)
        
        # 2006
        secondary_response = data[data['gene_names'].isin(['{Il12b}','{Il6}','{Nos2}','{Mx1}','{Mx2}','{Marco}','{Cmpk2}','{Rsad2}'])]
        delayed = data[data['gene_names'].isin(['{Ccl5}','{Saa3}','{Ifnb1}','{Ccl2}','{Ifit1}','{Ifit3}','{Peli1}','{Cxcl10}','{Traf1}'])]
        early_response = data[data['gene_names'].isin(['{Cxcl2}','{Tnf}','{Ptgs2}','{Egr1}','{Nfkbiz}','{Il1b}','{Tnfaip3}','{Cxcl1}','{Socs3}'])]
        
        # 2009
        swi_snf_ind_cpg = data[data['gene_names'].isin(['{Cd83}', '{Nr4a1}', '{Ccrn4l}', '{Irf1}', '{Nfkbiz}', '{Sod2}', '{Pim1}', '{Socs3}', '{Tnfaip3}', '{Egr1}', '{Junb}', '{Marcksl1}', '{Bcl3}', '{Icam1}', '{Zfp36}', '{Tnfsf9}', '{Nfkbia}', '{Fosb}', '{Egr2}', '{Fos}', '{Nfkbid}', '{Ptgs2}', '{Tlr2}', '{Cxcl1}', '{Cxcl2}', '{Tnf}'])]
        swi_snf_ind_non_cpg = data[data['gene_names'].isin(['{Ccrl2}', '{Traf1}', '{Cxcl11}', '{Clec4e}', '{Il1a}', '{Csf2}', '{Il23a}', '{Ccl3}', '{Gbp2}', '{Gbp1}', '{Il1b}'])]
        swi_snf_dep_non_irf3 = data[data['gene_names'].isin(['{Map3k8}', '{Serpine1}', '{Arhgef3}', '{Vcam1}', '{Saa3}', '{Ccl2}', '{Il10}', '{Irg1}', '{Ikbke}', '{Ccl12}', '{Mmp13}'])]
        swi_snf_dep_irf3 = data[data['gene_names'].isin(['{Peli1}', '{Ifit2}', '{Cxcl10}', '{Ifit1}', '{Ifnb1}', '{Ccl5}', '{Ifit3}',])]
        secondary = data[data['gene_names'].isin(['{Irf7}', '{Il6}', '{Il12b}', '{Nos2}', '{Lcn2}', '{Marco}', '{Mx1}', '{Mx2}', '{Serpinb3b}',])]
        
        # Cpg Islands
        refseq = data[data['gene_names'] != 0]
        high_cpg = refseq[refseq['has_cpg'] == 1]
        low_cpg = refseq[refseq['has_cpg'] == 0]
        
        supersets = ((secondary_response, delayed, early_response),
                     (swi_snf_ind_cpg, swi_snf_ind_non_cpg, swi_snf_dep_non_irf3, swi_snf_dep_irf3, secondary),
                     (high_cpg, low_cpg))
        superlabels = (('Secondary Response', 'Delayed Early\nResponse','Early Response'),
                       ('Early SWI/SNF Ind.\n CpG-Island Prom.', 'Early SWI/SNF Ind.\n No CpG-Island',
                        'Early SWI/SNF Dep.\n Irf3 Ind.', 'Early SWI/SNF Dep.\n Irf3 Dep.',
                        'Secondary\nSWI/SNF Dep.'),
                       ('Has CpG Island', 'No CpG Island'))
        title_groups = ('Ramirez-Carrozzi 2006', 'Ramirez-Carrozzi 2009', 'CpG Islands')
        for labels, sets, group in zip(superlabels, supersets, title_groups):
            for d in sets: print len(d)
            vals = [d['dex_over_kla_{0}_lfc'.format(rep)] for d in sets]
            
            title = 'Dex Transrepression for Genes Groups defined by {0}'.format(group)
            ax = yzer.boxplot(vals, labels, 
                             title=title, xlabel='Gene Group', 
                             ylabel='log2(KLA+Dex/KLA)', 
                             show_outliers=False, show_plot=False, 
                             save_dir=img_dirpath
                             )