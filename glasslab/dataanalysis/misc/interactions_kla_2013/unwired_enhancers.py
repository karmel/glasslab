'''
Created on Jan 2, 2013

@author: karmel

Many KLA interactions are just "rewired" from existing notx reactions.
Which are not?
'''

from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_sets')
    
    #img_dirpath = yzer.get_and_create_path(dirpath, 'boxplots_by_me2_in_notx')

    data = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_refseq_with_me2.txt'))
    data = data[data['count'] > 1]
    notx_enh = data[data['sequencing_run_id'] == 765]
    kla_30m_enh = data[data['sequencing_run_id'] == 766]
    kla_4h_enh = data[data['sequencing_run_id'] == 773]
    
    kla_30m_unwired = kla_30m_enh[~kla_30m_enh['id_2'].isin(notx_enh['id_2'])]
    print len(kla_30m_unwired)
    print kla_30m_unwired['ucsc_link_mm9']