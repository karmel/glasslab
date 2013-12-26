'''
Created on Dec 24, 2013

@author: karmel

Given a (Homer) peak file, merge all peaks that are within 12.5kb of each other,
as per Whyte 2012- "Master Transcription Factors and Mediator Establish 
Super-Enhancers at Key Cell Identity Genes."

Keep track of original peaks.

No attempt in the original paper seems to have been made to remove genic or 
promoter regions. Therefore, no attempt is made here either.

Output is saved to a file.

Note that we are assuming on the order of tens of thousands of peaks. Thus,
we will not waste must effort writing temp files, and will assume we can hold
all relevant data in memory. This may need to be changed if our data set size
changes.

'''
from __future__ import division
from pandas.io import parsers
import sys
from pandas.core.frame import DataFrame
import os
from collections import OrderedDict

class PeakMerger(object):
    distance = 12500
    
    def import_file(self, filename, separator='\t', header=True, 
                    index_col=None, skiprows=None):
        if header: header_row = 0
        else: header_row = None
        data = parsers.read_csv(filename, 
                                sep=separator, 
                                header=header_row, 
                                index_col=index_col,
                                skiprows=skiprows)
        
        return data
    
    def sort_by_chr_bp(self, data):
        '''
        Sort data by chromosome then basepair start position, ascending.
        '''
        data = data.sort(['chr','start'])
        return data
    
    def merge_peaks(self, data):
        '''
        For each chromosome, merge peaks within .distance basepairs
        of each other.
        '''
        current_chr = None
        current_row = None
        
        merged = []
        for _, row in data.iterrows():
            if current_chr != row['chr']:
                # Restart counting.
                current_chr = row['chr']
                current_row = self._extract_current_row(row)
            else:
                if row['start'] - current_row['end'] <= self.distance:
                    # These peaks should be merged.
                    current_row['end'] = row['end']
                    current_row['merged_ids'].append(row['#PeakID'])
                    current_row['tag_count'] += row['Normalized Tag Count'] 
                else:    
                    # Close off this row.
                    merged.append(current_row)
                    current_row = self._extract_current_row(row)
        
        merged = DataFrame(merged)        
        return merged
    
    def _extract_current_row(self, row):
        current_row = OrderedDict([('chr', row['chr']),
                       ('start', row['start']),
                       ('end', row['end']),
                       ('merged_ids', [row['#PeakID']]),
                       ('tag_count', row['Normalized Tag Count'])]
                       )
        return current_row
    
    def output_peaks(self, data, filename):
        output_file, _ = os.path.splitext(filename)
        output_file = output_file + '_merged.txt'
        
        data = data[['chr', 'start', 'end', 'tag_count', 'merged_ids']]
        data = data.sort(['tag_count'], ascending=False)
        
        data.to_csv(output_file, sep='\t',  
                    header=True, index=True)
        return data
    
if __name__ == '__main__':
    
    # Get filepath
    try: filename = sys.argv[1]
    except IndexError: filename = '/Users/karmel/GlassLab/Notes_and_Reports/Super-Enhancers/whyte_2012/th1_tbet_peaks.txt'
    
    merger = PeakMerger()
    
    data = merger.import_file(filename, skiprows=39)
    data = merger.sort_by_chr_bp(data)
    data = merger.merge_peaks(data)
    data = merger.output_peaks(data, filename)
    
    
    