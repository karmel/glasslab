'''
Created on Dec 8, 2012

@author: karmel

For loading in HiC interaction tags as filtered and paired by Homer.
'''
from __future__ import division
import os
import re
from glasslab.utils.scripting import GlassOptionParser
from optparse import make_option
import traceback

from django.db import connection, transaction
from glasslab.sequencing.datatypes.tag import GlassTag
from multiprocessing import Pool
from glasslab.config import current_settings
from glasslab.sequencing.pipeline.add_short_reads import check_input, _print,\
    create_schema
from glasslab.utils.database import execute_query_without_transaction

class FastqOptionParser(GlassOptionParser):
    options = [
               make_option('-f', '--file_name',action='store', type='string', dest='file_name', 
                           help='Path to Homer Tag directory for processing.'),
               make_option('-g', '--genome',action='store', type='string', dest='genome', default='mm9', 
                           help='Currently supported: mm8, mm8r, mm9, hg18, hg18r'),
               make_option('--project_name',action='store', type='string', dest='project_name',  
                           help='Optional name to be used as file prefix for created files.'),
               
               make_option('--schema_name',action='store', type='string', dest='schema_name',  
                           help='Optional name to be used as schema for created DB tables.'),
               
               make_option('--prep_table',action='store', dest='prep_table',
                           help='Skip transferring interactions from file to prep table; prep tag table will be used directly.'),
               
               ]
    
def copy_into_table_from_dir(tag_dir, f_name):
    if re.match(r'chr\d+\.tags\.tsv', f_name):
        _copy_into_table(tag_dir, f_name)
         
def _copy_into_table(tag_dir, f_name):
    full_path = os.path.join(tag_dir,f_name)
    f = file(full_path)
    try:
        connection.close()
        cursor = connection.cursor()
        cursor.copy_expert("""COPY "%s" (chromosome_1, "start_1",
                                        strand_1, count, length_1, chromosome_2, 
                                        "start_2", strand_2, length_2)
                                FROM STDIN WITH CSV DELIMITER E'\t'; """ % GlassTag.prep_table, f)
        transaction.commit_unless_managed()
    except Exception:
        _print('Encountered exception while trying to copy data:\n%s' % traceback.format_exc())
        raise
    
def upload_interaction_files(options, file_name, tag_dir):
    GlassTag.create_prep_table(file_name)
    
    file_names = os.listdir(tag_dir)
    processes = current_settings.ALLOWED_PROCESSES
    step_size = max(1,len(file_names)//processes)
    p = Pool(processes) 
    for start in xrange(0,len(file_names),step_size):
        try:
            p.apply_async(copy_into_table_from_dir,args=(tag_dir, file_names[start:(start + step_size)]))
        except Exception:
            raise Exception('Exception encountered while trying to upload tag files to tables. Traceback:\n%s'
                                % traceback.format_exc())
    p.close()
    p.join()
    
    connection.close()
    
def translate_prep_columns(file_name):
    '''
    Transfer prep tags to indexed, streamlined Glass tags for annotation.
    '''
    #GlassTag.set_table_name('tag_' + file_name)
    GlassTag.create_parent_table(file_name)
    GlassTag.create_partition_tables()
    GlassTag.translate_from_prep()
    GlassTag.add_record_of_tags()
        
def add_indices():
    # Execute after all the ends have been calculated,
    # as otherwise the insertion of ends takes far too long.
    GlassTag.add_indices()
    execute_query_without_transaction('VACUUM ANALYZE "%s";' % (GlassTag._meta.db_table))
        
if __name__ == '__main__':    
    run_from_command_line = True # Useful for debugging in Eclipse
    
    parser = FastqOptionParser()
    options, args = parser.parse_args()

    file_name = check_input(options)
     
    if not options.prep_table:
        _print('Creating schema if necessary.')
        create_schema()
        _print('Uploading tag file into table.')
        upload_interaction_files(options, file_name, options.file_name)
    else:
        _print('Skipping upload of interaction rows into table.')
        GlassTag.set_prep_table(options.prep_table)
    
    _print('Translating prep columns to integers.')
    translate_prep_columns(file_name)
    _print('Adding indices.')
    GlassTag.set_table_name('tag_' + file_name)
    add_indices()
