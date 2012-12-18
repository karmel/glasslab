#!/bin/bash
'''
Created on Sep 24, 2010

@author: karmel

A script capable of reading in a peak file and appropriately denormalizing into DB tables.

Run from the command line.

'''
from glasslab.utils.scripting import GlassOptionParser
from optparse import make_option
from glasslab.sequencing.datatypes.peak import GlassPeak
from glasslab.utils.parsing.delimited import DelimitedFileParser
from glasslab.sequencing.pipeline.add_short_reads import check_input, \
    create_schema, _print

class FastqOptionParser(GlassOptionParser):
    options = [
               make_option('-f', '--file_name',action='store', type='string', dest='file_name', 
                           help='Path to FASTQ file for processing.'),
               make_option('-o', '--output_dir',action='store', type='string', dest='output_dir'),
               make_option('-g', '--genome',action='store', type='string', dest='genome', default='mm9', 
                           help='Currently supported: mm8, mm8r, mm9, hg18, hg18r'),
               make_option('--project_name',action='store', type='string', dest='project_name',  
                           help='Optional name to be used as file prefix for created files.'),
            
               make_option('--schema_name',action='store', type='string', dest='schema_name',  
                           help='Optional name to be used as schema for created DB tables.'),
               
               make_option('--peak_type',action='store', dest='peak_type',  
                           help='What type of peak are we looking for? H4K3me1, PU_1, etc.? Should match a type in PeakType table.'),
               make_option('--not_homer',action='store', dest='not_homer', default=False, 
                           help='If the input file is not a HOMER peaks file, what is it? (bed, macs, sicer)'),
                           
               ]
    

def import_peaks(options, file_name, peaks_file_name, peak_type):
    '''
    Given a  peak file, save a temp table and store peak data.
    '''
    GlassPeak.create_table(file_name)
    #GlassPeak.set_table_name('peak_' + file_name)
    
    parser = DelimitedFileParser(peaks_file_name)
    parser.convert_line_endings()
    data = parser.get_array()
    if data[0][0] == 'chr':
        data = data[1:] # Header row
    
    for row in data:
        if not options.not_homer:
            peak = GlassPeak.init_from_homer_row(row)
        else:
            peak = getattr(GlassPeak, 'init_from_{0}_row'.format(options.not_homer))(row)

        peak.save()

    if options.not_homer and peak_type.diffuse:
        GlassPeak.score_sicer_peaks()
    
    GlassPeak.add_indices()
    GlassPeak.add_record_of_tags(peak_type=peak_type)
    
    
if __name__ == '__main__': 
    run_from_command_line = True # Useful for debugging in Eclipse
    
    parser = FastqOptionParser()
    options, args = parser.parse_args()
    
    # Allow for easy running from Eclipse
    if not run_from_command_line:
        options.file_name = '/Users/karmel/GlassLab/SourceData/ThioMac_Lazar/test'
        options.output_dir = 'first_test'
        options.peak_table = 'enriched_peaks_enriched_peaks_first_test_2010-09-30_16-04-52_454502'
    
    file_name = check_input(options)
    
    peaks_file_name = options.file_name
    
    peak_type = GlassPeak.peak_type(options.peak_type or options.project_name)
    _print('Creating schema if necessary.')
    create_schema()
    _print('Saving peaks to table.')
    import_peaks(options, file_name, peaks_file_name, peak_type)
