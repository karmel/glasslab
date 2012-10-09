'''
Created on Nov 8, 2010

@author: karmel
'''
from glasslab.glassatlas.datatypes.transcript import CellTypeBase
from glasslab.utils.scripting import GlassOptionParser
from optparse import make_option
from glasslab.config import current_settings
from glasslab.sequencing.datatypes.peak import GlassPeak

class TranscriptsFromTagsParser(GlassOptionParser):
    options = [
               make_option('-c', '--cell_type',action='store', type='string', dest='cell_type', 
                           help='Cell type for this run? Options are: %s' % ','.join(CellTypeBase.get_correlations().keys())),
               make_option('-t', '--peak_table',action='store', type='string', dest='peak_table', 
                           help='Table name from which to load peaks. Appended to schema if schema is included. Otherwise used as is.'),
               make_option('-s', '--schema_name',action='store', type='string', dest='schema_name',  
                           help='Optional name to be used as schema for created DB tables.'),
               
                ]
if __name__ == '__main__':
    parser = TranscriptsFromTagsParser()
    options, args = parser.parse_args()
    
    run_from_cammand_line = True 
    if not run_from_cammand_line:
        options.schema_name = 'thiomac_groseq_nathan_2010_10'
        options.peak_table = 'peak_ncor_ko_kla_1h'
    
    if options.cell_type: current_settings.CELL_TYPE = options.cell_type
    cell_base = CellTypeBase().get_cell_type_base(current_settings.CELL_TYPE)()
    
    GlassPeak._meta.db_table = options.schema_name and '%s"."%s' % (options.schema_name, options.peak_table) \
                                or options.peak_table
    cell_base.peak_feature.add_from_peaks(GlassPeak._meta.db_table)

    
