'''
Created on Nov 8, 2010

@author: karmel
'''
from glasslab.glassatlas.datatypes.transcript import CellTypeBase
from glasslab.sequencing.datatypes.tag import GlassTag
from glasslab.utils.scripting import GlassOptionParser
from optparse import make_option
from glasslab.config import current_settings
from glasslab.utils.database import restart_server, discard_temp_tables

class TranscriptsFromTagsParser(GlassOptionParser):
    options = [
               make_option('-g', '--genome',action='store', type='string', dest='genome', default='mm9', 
                           help='Currently supported: mm8, mm8r, mm9, hg18, hg18r'),
               make_option('-c', '--cell_type',action='store', type='string', dest='cell_type', 
                           help='Cell type for this run? Options are: %s' % ','.join(CellTypeBase.get_correlations().keys())),
               make_option('-t', '--tag_table',action='store', type='string', dest='tag_table', 
                           help='Table name from which to load tags. Appended to schema if schema is included. Otherwise used as is.'),
               make_option('-s', '--schema_name',action='store', type='string', dest='schema_name',  
                           help='Optional name to be used as schema for created DB tables.'),
               make_option('-o', '--output_dir',action='store', type='string', dest='output_dir',  
                           help='Output directory for bed file.'),
               make_option('-p','--processes',action='store', dest='processes', default=None,  
                           help='How many processes can be used?'),
               make_option('--stitch_processes',action='store', dest='stitch_processes', default=None,  
                           help='How many processes can be used during stitching, which is memory-intensive?'),
               
               make_option('--stitch',action='store_true', dest='stitch',  
                           help='Should transcripts be stitched together?'),
               make_option('--set_density',action='store_true', dest='set_density',  
                           help='Set density? Performed with tag stitching, or separately if --stitch=False.'),
               make_option('--draw_edges',action='store_true', dest='draw_edges',  
                           help='Should the edges between transcripts be created and saved?'),
               make_option('--score',action='store_true', dest='score',  
                           help='Should transcripts be stored?'),
               make_option('--set_start_end_tss',action='store_true', dest='set_start_end_tss',  
                           help='Should transcripts be updated with either whole start_end values or just the TSS in the field set_start_end_tss?'),
               make_option('--no_extended_gaps',action='store_true', dest='no_extended_gaps',  
                           help='Should extended gaps (i.e., under RefSeq regions) be allowed?'),
               make_option('--staging',action='store_true', dest='staging', default=False,  
                           help='Use the transcript database with the suffix _staging?'),
               make_option('--max_edge',action='store', dest='max_edge',  
                           help='What value of MAX_EDGE should be used for drawing transcript edges?'),
                ]

if __name__ == '__main__':
    parser = TranscriptsFromTagsParser()
    options, args = parser.parse_args()
    
    run_from_cammand_line = True 
    if not run_from_cammand_line:
        options.schema_name = 'thiomac_groseq_nathan_2010_10'
        options.tag_table = 'tag_ncor_ko_kla_1h'
    
    if options.processes:
        current_settings.ALLOWED_PROCESSES = int(options.processes)
    
    if options.cell_type: current_settings.CELL_TYPE = options.cell_type
    cell_base = CellTypeBase().get_cell_type_base(current_settings.CELL_TYPE)()
    
    allow_extended_gaps = True
    if options.no_extended_gaps: allow_extended_gaps = False
    
    if options.staging: current_settings.STAGING = current_settings.STAGING_SUFFIX

    if options.tag_table:
        GlassTag._meta.db_table = options.schema_name and '%s"."%s' % (options.schema_name, options.tag_table) \
                                    or options.tag_table
        cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        #cell_base.glass_transcript.force_vacuum_prep()
        #print 'Restarting server...'
        #restart_server()
    
    if options.stitch_processes:
        curr_processes = current_settings.ALLOWED_PROCESSES 
        current_settings.ALLOWED_PROCESSES = int(options.stitch_processes)
    if options.stitch:
        cell_base.glass_transcript.stitch_together_transcripts(
                        allow_extended_gaps=allow_extended_gaps, 
                        extension_percent=options.extension_percent, 
                        set_density=options.set_density)
        restart_server()
    elif options.set_density:
        cell_base.glass_transcript.set_density(allow_extended_gaps=allow_extended_gaps,
                                               extension_percent=options.extension_percent)

    if options.stitch_processes:
        current_settings.ALLOWED_PROCESSES = curr_processes
   
    if options.draw_edges:
        cell_base.glass_transcript.draw_transcript_edges()
    
    if options.score:
        cell_base.glass_transcript.set_scores()
        
    if options.set_start_end_tss:
        cell_base.glass_transcript.set_start_end_tss()
        
    if options.output_dir:
        cell_base.filtered_glass_transcript.generate_bed_file(options.output_dir)
        
    discard_temp_tables()
    
