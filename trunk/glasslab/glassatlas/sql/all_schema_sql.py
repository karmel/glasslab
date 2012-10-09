'''
Created on Feb 23, 2011

@author: karmel
'''
import sys
from glasslab.glassatlas.sql.transcripts_from_tags_functions import sql as transcript_function_sql
from glasslab.glassatlas.sql.transcripts_from_prep_functions import sql as transcript_from_prep_function_sql
from glasslab.glassatlas.sql.features_functions import sql as features_function_sql
from glasslab.glassatlas.sql.glassatlas_table_generator import GlassAtlasTableGenerator


if __name__ == '__main__':
    genome = len(sys.argv) > 1 and sys.argv[1] or 'gap3_100_10'
    cell_type = len(sys.argv) > 2 and sys.argv[2] or 'thiomac'
    subset = len(sys.argv) > 3 and sys.argv[3] or False
    
    generator = GlassAtlasTableGenerator(genome=genome, cell_type=cell_type)
    if subset == 'final':
        print generator.final_set()
        print transcript_from_prep_function_sql(genome, cell_type)
        print features_function_sql(genome, cell_type)
    else: 
        print generator.all_sql()
        print transcript_function_sql(genome, cell_type)
        print transcript_from_prep_function_sql(genome, cell_type)
        print features_function_sql(genome, cell_type)
