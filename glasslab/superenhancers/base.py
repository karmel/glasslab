'''
Created on Dec 27, 2013

@author: karmel
'''
from pandas.io import parsers


class PeakFinder(object):
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