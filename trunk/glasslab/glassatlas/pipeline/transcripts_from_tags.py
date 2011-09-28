'''
Created on Nov 8, 2010

@author: karmel
'''
from glasslab.glassatlas.datatypes.transcript import CellTypeBase
from glasslab.sequencing.datatypes.tag import GlassTag
from glasslab.utils.scripting import GlassOptionParser
from optparse import make_option
from glasslab.config import current_settings

class TranscriptsFromTagsParser(GlassOptionParser):
    options = [
               make_option('-c', '--cell_type',action='store', type='string', dest='cell_type', 
                           help='Cell type for this run? Options are: %s' % ','.join(CellTypeBase.get_correlations().keys())),
               make_option('-t', '--tag_table',action='store', type='string', dest='tag_table', 
                           help='Table name from which to load tags. Appended to schema if schema is included. Otherwise used as is.'),
               make_option('-s', '--schema_name',action='store', type='string', dest='schema_name',  
                           help='Optional name to be used as schema for created DB tables.'),
               make_option('-o', '--output_dir',action='store', type='string', dest='output_dir',  
                           help='Output directory for bed file.'),
               make_option('--remove_rogue_run',action='store_true', dest='remove_rogue_run',  
                           help='Should the records from this run be removed?'),
               make_option('--skip_stitching',action='store_true', dest='skip_stitching',  
                           help='Should the stitching together of transcripts be skipped?'),
               make_option('--set_density',action='store_true', dest='set_density',  
                           help='Should the stitching together of transcripts be replaced with forcing reset of average tags?'),
               make_option('--draw_edges',action='store_true', dest='draw_edges',  
                           help='Should the edges between transcripts be created and saved?'),
               make_option('--reset_score_thresholds',action='store_true', dest='reset_score_thresholds',  
                           help='Should the expected score thresholds for each chromosome and strand be reset?'),
               make_option('--skip_scoring',action='store_true', dest='skip_scoring',  
                           help='Should the scoring of transcripts be skipped?'),
               make_option('--associate_nucleotides',action='store_true', dest='associate_nucleotides',  
                           help='Should obtaining nucleotide sequences to transcripts be attempted?'),
               make_option('--reassociate',action='store_true', dest='reassociate',  
                           help='Should reassociation of peak features to transcripts be attempted?'),
               make_option('--no_extended_gaps',action='store_true', dest='no_extended_gaps',  
                           help='Should extended gaps (i.e., under RefSeq regions) be allowed?'),
               make_option('--staging',action='store_true', dest='staging', default=False,  
                           help='Use the transcript database with the suffix _staging?'),
                ]

if __name__ == '__main__':
    parser = TranscriptsFromTagsParser()
    options, args = parser.parse_args()
    
    run_from_cammand_line = True 
    if not run_from_cammand_line:
        options.schema_name = 'thiomac_groseq_nathan_2010_10'
        options.tag_table = 'tag_ncor_ko_kla_1h'
    
    if options.cell_type: current_settings.CURRENT_CELL_TYPE = options.cell_type
    cell_base = CellTypeBase().get_cell_type_base(current_settings.CURRENT_CELL_TYPE)()
    
    allow_extended_gaps = True
    if options.no_extended_gaps: allow_extended_gaps = False
    
    if options.staging: current_settings.STAGING = current_settings.STAGING_SUFFIX

    #cell_base.glass_transcript.turn_off_autovacuum()    
    if options.tag_table:
        GlassTag._meta.db_table = options.schema_name and '%s"."%s' % (options.schema_name, options.tag_table) \
                                    or options.tag_table
        cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        #cell_base.glass_transcript.force_vacuum_prep()
    elif options.remove_rogue_run:
        cell_base.glass_transcript.remove_rogue_run()
        cell_base.glass_transcript.force_vacuum_prep()
    
    if not options.skip_stitching:
        cell_base.glass_transcript.stitch_together_transcripts(
                        allow_extended_gaps=allow_extended_gaps, set_density=options.set_density)
        cell_base.glass_transcript.force_vacuum_prep()
    elif options.set_density:
        cell_base.glass_transcript.set_density(allow_extended_gaps=allow_extended_gaps)
        cell_base.glass_transcript.force_vacuum_prep()
    
    if options.draw_edges:
        cell_base.glass_transcript.draw_transcript_edges()
        #cell_base.glass_transcript.force_vacuum()
    
    if options.associate_nucleotides:
        cell_base.glass_transcript.associate_nucleotides()
    
    if options.reset_score_thresholds:
        cell_base.glass_transcript.set_score_thresholds()
        
    if not options.skip_scoring:
        cell_base.glass_transcript.set_scores()
        
    if options.reassociate:
        cell_base.peak_feature.update_peak_features_by_transcript()
        cell_base.glass_transcript.mark_all_reloaded()
        
    if options.output_dir:
        cell_base.filtered_glass_transcript.generate_bed_file(options.output_dir)
        
    #cell_base.glass_transcript.turn_on_autovacuum()
    
