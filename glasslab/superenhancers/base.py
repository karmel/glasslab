'''
Created on Dec 27, 2013

@author: karmel
'''
from pandas.io import parsers


class PeakFinder(object):
    def find_header(self, filename):
        '''
        For a HOMER file, find the line number of the header.
        '''
        
        for i, line in enumerate(file(filename).readlines()):
            if line[:5] == '#Peak': return i
        return 0
    
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