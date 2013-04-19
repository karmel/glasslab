'''
Created on Apr 12, 2013

@author: karmel
'''
from __future__ import division
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer


if __name__ == '__main__':
    yzer = MotifAnalyzer()
    
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/CD4TCells/TReg_enhancers/2013_04_01'
    dirpath = yzer.get_path(dirpath)
    data = yzer.import_file(yzer.get_filename(dirpath, 'th1_with_stat1_ko.txt')).fillna(0)
    
    print len(data)
    data = data[data['tss_id'] == 0]
    data['ko_ratio'] = data['ko_id']/data['th1_id']
    data['treg_ratio'] = data['treg_id']/data['th1_id']
    enh = len(data)
    print enh
    print sum(data['ko_ratio'] < .5), sum(data['ko_ratio'] < .5)/enh
    print sum(data['treg_ratio'] < .5), sum(data['treg_ratio'] < .5)/enh
    print sum((data['treg_ratio'] < .5) & (data['ko_ratio'] < .5))
    print sum((data['treg_ratio'] < .5) & (data['ko_ratio'] < .5))/enh
                    