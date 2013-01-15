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

    data = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_enhancer_with_anything.txt'))
    data = data[data['count'] > 1]
    notx_enh = data[data['sequencing_run_id'] == 765]
    kla_30m_enh = data[data['sequencing_run_id'] == 766]
    kla_4h_enh = data[data['sequencing_run_id'] == 773]
    
    print len(kla_4h_enh['id_2'].unique())
    kla_4h_unwired = kla_4h_enh[~kla_4h_enh['id_2'].isin(notx_enh['id_2'])]
    print len(kla_4h_unwired['id_2'].unique())
    print kla_4h_unwired['ucsc_link_mm9']