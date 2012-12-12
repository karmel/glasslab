'''
Created on Nov 8, 2010

@author: karmel
'''
from glasslab.glassatlas.datatypes.transcript import CellTypeBase
from glasslab.utils.scripting import GlassOptionParser
from optparse import make_option
from glasslab.config import current_settings
from glasslab.utils.database import discard_temp_tables

class InteractionsParser(GlassOptionParser):
    options = [
               make_option('-g', '--genome',action='store', type='string', dest='genome', default='mm9', 
                           help='Currently supported: mm8, mm8r, mm9, hg18, hg18r'),
               make_option('-c', '--cell_type',action='store', type='string', dest='cell_type', 
                           help='Cell type for this run? Options are: %s' % ','.join(CellTypeBase.get_correlations().keys())),
               make_option('-t', '--interaction_table',action='store', type='string', dest='interaction_table', 
                           help='Table name from which to load interactions. Appended to schema if schema is included. Otherwise used as is.'),
               make_option('-s', '--schema_name',action='store', type='string', dest='schema_name',  
                           help='Optional name to be used as schema for created DB tables.'),
               
                ]
if __name__ == '__main__':
    parser = InteractionsParser()
    options, args = parser.parse_args()
    
    if options.cell_type: current_settings.CELL_TYPE = options.cell_type
    cell_base = CellTypeBase().get_cell_type_base(current_settings.CELL_TYPE)()
    
    source_table = options.schema_name and '{0}"."{1}'.format(options.schema_name, options.interaction_table) \
                                or options.interaction_table
    cell_base.glass_transcript.associate_interactions(source_table)

    discard_temp_tables();
