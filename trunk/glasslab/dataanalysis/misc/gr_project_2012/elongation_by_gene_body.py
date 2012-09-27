'''
Created on Sep 26, 2012

@author: karmel
'''

from glasslab.dataanalysis.misc.gr_project_2012.width_buckets import get_data_with_bucket_score
from glasslab.dataanalysis.base.datatypes import TranscriptAnalyzer

if __name__ == '__main__':
    learner = TranscriptAnalyzer()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/classification'
    dirpath = learner.get_path(dirpath)
    
    
    # First time file setup
    data = get_data_with_bucket_score(learner, dirpath)
    
    print data.columns