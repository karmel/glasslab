'''
Created on Mar 22, 2012

@author: karmel
'''
from __future__ import division
from pandas.io import parsers
from pandas.core.frame import DataFrame
import os

class TranscriptAnalyzer(object):
    '''
    Using pandas, does some basic importing and filtering of transcript vectors.
    '''    
    def get_path(self, dirpath=''):
        locs = ('','/Users','/Volumes')
        for loc in locs:
            path = os.path.join(loc,dirpath)
            if os.path.exists(path): return path
        raise IOError('Could not find file %s.' % dirpath)
    
    def get_filename(self, *args):
        return os.path.join(*args)
    
    def get_and_create_path(self, *args):
        dirpath = os.path.join(*args)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        return dirpath
    
    def import_file(self, filename, separator='\t', header=True, index_col=None, skiprows=None):
        if header: header_row = 0
        else: header_row = None
        data = parsers.read_csv(filename, 
                                sep=separator, 
                                header=header_row, 
                                index_col=index_col,
                                skiprows=skiprows)
        
        return data
    
    def get_refseq(self, data, colname='has_refseq'):
        return data[data[colname] != 0]
    
    def normalize(self, data, colname, normfactor, suffix='_norm'):
        data[colname + suffix] = data[colname]*normfactor
        return data
    
    def get_gene_list(self, data, colname='gene_names'):
        return map(lambda x: x.replace('{','').replace('}',''), data[colname].values)
    
    def get_gene_names(self, data, colname='gene_names', add_quotes=False):
        joiner = add_quotes and "', '" or ', '
        return joiner.join(data[colname].values).replace('{','').replace('}','')

    def collapse_strands(self, data):
        '''
        If the data is strand specific, we may want to collapse
        overlapping sense and antisense transcripts for motif finding.
        
        This is relevant in the case of eRNA, where transcription is 
        bidirectional, and we want the center of transcription.
        
        Testing:
        data = DataFrame({'chr_name': ['chr1', 'chr1', 'chr1','chr2','chr2',], 
                       'transcription_start': [0, 10, 21, 0, 5],
                       'transcription_end': [10, 20, 30, 10, 7], 
                       'strand': [0, 1, 0, 1, 0]},
                     )
        should yield
              chr_name strand transcription_end transcription_start
                0     chr1      0                20                   0
                2     chr1      0                30                  21
                3     chr2      1                10                   0

        '''
        compressed = []
        try: 
            chr_key = 'chr_name'
            ordered = data.sort_index(by=[chr_key,'transcription_start'])
        except KeyError: 
            chr_key = 'chromosome_id'
            ordered = data.sort_index(by=[chr_key,'transcription_start'])
        last = None
        for _, trans in ordered.iterrows():
            try:
                if trans[chr_key] == last[chr_key]\
                    and trans['transcription_start'] <= last['transcription_end'] + 200:
                        last['transcription_end'] = max(trans['transcription_end'],last['transcription_end'])
                        for key in last.index:
                            if key.find('_start') >= 0:
                                last[key] = min(last[key], trans[key])
                            elif key.find('_end') >= 0:
                                last[key] = max(last[key], trans[key])
                            elif key.find('id') == 0 or key.find('_id') >= 0:
                                pass # Keep single id
                            else: 
                                # Average any other vals, list any other strings
                                try: last[key] = (last[key] + trans[key])/2
                                except TypeError: last[key] = ','.join([last[key], trans[key]])
                else: 
                    compressed.append((last.name, last))
                    last = trans
            except TypeError: last = trans
        compressed.append((last.name, last))
            
        return DataFrame(dict(compressed)).transpose()