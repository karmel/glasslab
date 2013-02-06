'''
Created on Jan 29, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/Glass Atlas/Demo-data'
    dirpath = yzer.get_path(dirpath)
    
    #img_dirpath = yzer.get_and_create_path(dirpath, 'distal_remainder')
    
    data = yzer.import_file(yzer.get_filename(dirpath,'distal_transcript_remainder.txt'))
    data['length'] = data['transcription_end'] - data['transcription_start'] + 1
    data = data.sort('length', ascending=False)
    print data.columns.values
    for i, row in data[:10].iterrows(): print row.values

    print data['length'].describe()
    
    near_refseq = yzer.import_file(yzer.get_filename(dirpath,'distal_transcript_remainder_near_refseq.txt'))
    near_me2 = yzer.import_file(yzer.get_filename(dirpath,'me2_peaks_within_2000bp.txt'))
    near_refseq_not_me2 = near_refseq[~near_refseq['id'].isin(near_me2['id'])]
    print len(near_refseq_not_me2)