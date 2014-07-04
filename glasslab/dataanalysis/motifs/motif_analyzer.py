'''
Created on Mar 22, 2012

@author: karmel
'''
from glasslab.dataanalysis.base.datatypes import TranscriptAnalyzer
import os
import subprocess
import math

class MotifAnalyzer(TranscriptAnalyzer):
    '''
    For use with the HOMER motif finding program for analysis
    of enriched motifs in transcript sets.
    '''
    homer_bin = '/Applications/bioinformatics/homer/bin/'
    default_length = [8, 10, 12] # For use with findMotifsGenome.pl
    default_size = 200
    
    
    def run_homer(self, data, project_name, dirpath,
                  center=True, reverse=False, preceding=False,
                  size=None, length=None, number_of_motifs=20,
                  bg='', genome='mm9', cpus=1,
                  files_already_prepped=False, mock=False):
        '''
        Run HOMER on passed set of transcripts, assuming a transcription_start,
        transcription_end, and chr_name.
        '''
        fullpath, homer_filename = self.prep_files_for_homer(data,
                            project_name, dirpath, 
                            center=center, reverse=reverse, preceding=preceding, size=size,
                            files_already_prepped=files_already_prepped)
        
        self.run_homer_with_pos_file(homer_filename, fullpath, center, size,
                                     length, number_of_motifs, bg, genome, cpus, 
                                     mock=mock)
        
    
    def prep_files_for_homer(self, data_orig, project_name, dirpath,
                  center=False, reverse=False, preceding=False, size=None,
                  files_already_prepped=False):
        '''
        Create tab-delimited files with Windows carriage returns for Homer.
        '''
        fullpath = self.get_filename(dirpath, project_name)
        self.make_directory(fullpath, check_exists=False)
        
        # Create HOMER-compatible file.
        region_filename = self.get_filename(fullpath, project_name + '_regions.txt')
        homer_filename = self.get_filename(fullpath, project_name + '_regions_for_homer.txt')
        
        if not files_already_prepped:
            # We don't want to modify the data.
            data = data_orig.copy()
            
            if not size: size = self.default_size
            
            try: data['strand'] = data['strand'].apply(int)
            except KeyError: data['strand'] = 0
            
            # Which key should we use?
            start_key = 'transcription_start'
            end_key = 'transcription_end'
            try: data[start_key]
            except KeyError: start_key = 'start'
            try: data[end_key]
            except KeyError: end_key = 'end'

            if not center:
                first_strand, second_strand = 0, 1
                if reverse:
                    # Take the end of the transcript, rather than the beginning. 
                    first_strand, second_strand = 1, 0
                if preceding:
                    # First, expand the whole transcript on either side
                    data[start_key] = data[start_key] - size
                    data[end_key] = data[end_key] + size
                    # Then grab strand-dependent beginning and end
                    data['transcription_end_alt'] = data[data['strand'] == first_strand][start_key] + size 
                    data['transcription_start_alt'] = data[data['strand'] == second_strand][end_key] - size 
                else:
                    data['transcription_end_alt'] = data[data['strand'] == first_strand][start_key] + size 
                    data['transcription_start_alt'] = data[data['strand'] == second_strand][end_key] - size
                data[start_key] = data['transcription_start_alt'].fillna(data[start_key])
                data[end_key] = data['transcription_end_alt'].fillna(data[end_key])
            data[start_key] = data[start_key].apply(int)
            data[end_key] = data[end_key].apply(int)
            
            
            # Pandas keeps pulling wrong columns? Subsetting data explicitly.
            data_subset = data[['id', 'chr_name', start_key, end_key, 'strand']]
            data_subset.to_csv(region_filename,
                        header=False, index=False, sep='\t')
            # Convert to Windows line endings
            subprocess.check_call("sed ':a;N;$!ba;s/\n/\r/g' %s > %s" % (
                                                        self.sanitize_filename(region_filename),
                                                        self.sanitize_filename(homer_filename)),
                                                        shell=True)
            
        return fullpath, homer_filename
        
    def run_homer_with_pos_file(self, homer_filename, dirpath,
                                center=True, size=None, length=None, number_of_motifs=10,
                                bg='', genome='mm9', cpus=1, mock=False):
        
        
        if not size: size = self.default_size
        if center: size = '-%d,+%d' % (int(math.floor(.5 * size)), int(math.ceil(.5 * size)))
        else: size = str(size)
        
        if not length: length = self.default_length
        length = ','.join(map(str, length))
        
        if bg: bg = '-bg "%s"' % bg
        else: bg = ''
        
        fullpath = self.get_filename(dirpath, 'homer_motifs_size_%s_len_%s' % (
                                        size.replace(',', '_'), length.replace(',', '-'))
                                )
        if not mock: self.make_directory(fullpath)
        
        command = '%sfindMotifsGenome.pl %s %s "%s" -size %s -len %s -p %d -S %d %s' % (
                                        self.homer_bin, self.sanitize_filename(homer_filename),
                                        genome, fullpath,
                                        size, length, cpus, number_of_motifs, bg)
        if not mock: subprocess.check_call(command, shell=True)
        
        print 'Successfully executed command:\n%s' % command
        
    def make_directory(self, fullpath, check_exists=True):
        if check_exists: self.check_exists(fullpath)
        
        try: os.mkdir(fullpath)
        except OSError, e: 
            if check_exists: raise e
    
    def check_exists(self, fullpath):
        if os.path.exists(fullpath): raise Exception('Directory {0} already exists!'.format(fullpath))
        
    def sanitize_filename(self, filename):
        return filename.replace(' ', '\ ')
    

    
    
