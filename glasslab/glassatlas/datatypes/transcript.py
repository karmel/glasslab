'''
Created on Nov 8, 2010

@author: karmel
'''
from __future__ import division
from glasslab.config import current_settings
from django.db import models, connection, utils
from glasslab.genomereference.datatypes import Chromosome,\
    SequenceTranscriptionRegion, InfrastructureTranscriptionRegion,\
    NonCodingTranscriptionRegion, SequencingRun
from glasslab.utils.datatypes.basic_model import BoxField, GlassModel
from multiprocessing import Pool
from glasslab.utils.database import execute_query,\
    execute_query_without_transaction, fetch_rows, discard_temp_tables
import os
from random import randint
from django.db.models.aggregates import Max
from datetime import datetime
import traceback

# The tags returned from the sequencing run are shorter than we know them to be biologically
# We can therefore extend the mapped tag region by a set number of bp if an extension is passed in
TAG_EXTENSION = 0

MAX_GAP = 0 # Max gap between transcripts from the same run
MAX_STITCHING_GAP = MAX_GAP # Max gap between transcripts being stitched together
MAX_EDGE = 20 # Max edge length of transcript graph subgraphs to be created
EDGE_SCALING_FACTOR = 20 # Number of transcripts per DENSITY_MULTIPLIER bp required to get full allowed edge length
DENSITY_MULTIPLIER = 1000 # Scaling factor on density-- think of as bps worth of tags to consider
MIN_SCORE = 6 # Hide transcripts with scores below this threshold.

def multiprocess_all_chromosomes(func, cls, *args, **kwargs):
    ''' 
    Convenience method for splitting up queries based on glass tag id.
    '''
    processes = current_settings.ALLOWED_PROCESSES
    
    if not current_settings.CHR_LISTS:
        try:
            try:
                # Note that we accept a kwarg use_table
                all_chr = fetch_rows('''
                    SELECT chromosome_id as id
                    FROM "{0}" 
                    GROUP BY chromosome_id ORDER BY COUNT(chromosome_id) DESC;'''.format(
                                    kwargs.get('use_table',None) or cls._meta.db_table))
            except utils.DatabaseError:
                # Prep table instead?
                all_chr = fetch_rows('''
                    SELECT chromosome_id as id
                    FROM "{0}" 
                    GROUP BY chromosome_id ORDER BY COUNT(chromosome_id) DESC;'''.format(
                                    getattr(cls,'prep_table',None)
                                        or cls.cell_base.glass_transcript_prep._meta.db_table))
                
            all_chr = zip(*all_chr)[0]
            if not all_chr: raise Exception
            
        except Exception:
            # cls in question does not have explicit relation to chromosomes; get all
            all_chr = current_settings.GENOME_CHOICES[current_settings.GENOME]['chromosomes']

        # Chromosomes are sorted by count descending, so we want to snake them
        # back and forth to create even-ish groups. 
        chr_sets = [[] for _ in xrange(0, processes)]
        for i,chrom in enumerate(all_chr):
            if i and not i % processes: chr_sets.reverse()
            chr_sets[i % processes].append(chrom)
        
        # Reverse every other group to even out memory requirements.
        for i, chr_set in enumerate(chr_sets):
            if i % 2 == 0: chr_set.reverse()
            
        current_settings.CHR_LISTS = chr_sets
        print 'Determined chromosome sets:\n{0}'.format(str(current_settings.CHR_LISTS))
    
    p = Pool(processes)
    for chr_list in current_settings.CHR_LISTS:
        p.apply_async(func, args=[cls, chr_list,] + list(args))
    p.close()
    p.join()

# The following methods wrap bound methods. This is necessary
# for use with multiprocessing. Note that getattr with dynamic function names
# doesn't seem to work either.
def wrap_errors(func, *args):
    try: func(*args)
    except Exception:
        print 'Encountered exception in wrapped function:\n{0}'.format(traceback.format_exc())
        raise
   
def wrap_add_transcripts_from_groseq(cls, chr_list, *args): wrap_errors(cls._add_transcripts_from_groseq, chr_list, *args)
def wrap_stitch_together_transcripts(cls, chr_list, *args): wrap_errors(cls._stitch_together_transcripts, chr_list, *args)
def wrap_set_density(cls, chr_list, *args): wrap_errors(cls._set_density, chr_list, *args)
def wrap_draw_transcript_edges(cls, chr_list): wrap_errors(cls._draw_transcript_edges, chr_list)
def wrap_set_scores(cls, chr_list): wrap_errors(cls._set_scores, chr_list)
def wrap_set_start_end_tss(cls, chr_list): wrap_errors(cls._set_start_end_tss, chr_list)
def wrap_associate_interactions(cls, chr_list, *args): wrap_errors(cls._associate_interactions, chr_list, *args)
def wrap_force_vacuum(cls, chr_list): wrap_errors(cls._force_vacuum, chr_list)

class CellTypeBase(object):
    cell_type = ''
    
    @classmethod
    def get_correlations(cls):
        from glasslab.glassatlas.datatypes.celltypes.thiomac import ThioMacBase
        from glasslab.glassatlas.datatypes.celltypes.bmdc import BMDCBase
        from glasslab.glassatlas.datatypes.celltypes.cd4tcell import CD4TCellBase
        from glasslab.glassatlas.datatypes.celltypes.default import DefaultBase, RefSeqBase
        return {'default': DefaultBase,
                'thiomac': ThioMacBase,
                'cd4tcell': CD4TCellBase,
				'bmdc':BMDCBase,
                'refseq': RefSeqBase}
      
    @property
    def glass_transcript(self): return GlassTranscript
    @property
    def glass_transcript_prep(self): return GlassTranscriptPrep
    @property
    def filtered_glass_transcript(self): return FilteredGlassTranscript
    @property
    def glass_transcript_source(self): return GlassTranscriptSource
    @property
    def glass_transcript_source_prep(self): return GlassTranscriptSourcePrep
    @property
    def glass_transcript_sequence(self): return GlassTranscriptSequence
    @property
    def glass_transcript_non_coding(self): return GlassTranscriptNonCoding
    @property
    def peak_feature(self): 
        from glasslab.glassatlas.datatypes.feature import PeakFeature
        return PeakFeature
    
    def get_transcript_models(self):
        return [self.glass_transcript, self.glass_transcript_prep, 
                self.glass_transcript_source, self.glass_transcript_source_prep, 
                self.glass_transcript_sequence, self.glass_transcript_non_coding,
                self.peak_feature]

    def get_cell_type_base(self, cell_type):
        correlations = self.__class__.get_correlations()
        try: return correlations[cell_type.lower()]
        except KeyError:
            raise Exception('Could not find models to match cell type {0}.'.format(cell_type)
                            + '\nOptions are: {0}'.format(','.join(correlations.keys())))
            
class TranscriptModelBase(GlassModel):
    cell_base = CellTypeBase()
    schema_base = 'glass_atlas_{0}_{1}'
    class Meta:
        abstract = True
        
class TranscriptionRegionBase(TranscriptModelBase):
    # Use JS for browser link to auto-include in Django Admin form 
    ucsc_session_args   = '''hgS_doOtherUser=submit&amp;hgS_otherUserName=Karmel&amp;hgS_otherUserSessionName=%s''' % current_settings.UCSC_SESSION
    ucsc_browser_link_1 = '''<a href="#" onclick="window.open('http://genome.ucsc.edu/cgi-bin/hgTracks?'''
    ucsc_browser_link_2 = '''&amp;db=''' + current_settings.GENOME + '''&amp;position=' + '''\
                        + ''' django.jQuery('#id_chromosome').attr('title') '''\
                        + ''' + '%3A+' + django.jQuery('#id_transcription_start').val() '''\
                        + ''' + '-' + django.jQuery('#id_transcription_end').val(),'Glass Atlas UCSC View ' + '''\
                        + str(randint(100,99999)) + '''); return false;">'''
                        
    ucsc_browser_link = ucsc_browser_link_1 + ucsc_session_args + ucsc_browser_link_2\
                        + '''View in UCSC Browser</a>'''
                        
    chromosome              = models.ForeignKey(Chromosome, help_text=ucsc_browser_link)
    strand                  = models.IntegerField(max_length=1, help_text='0 for +, 1 for -')
    transcription_start     = models.IntegerField(max_length=12)
    transcription_end       = models.IntegerField(max_length=12, help_text='<span id="length-message"></span>')
    
    start_end               = BoxField(null=True, default=None, help_text='This is a placeholder for the PostgreSQL box type.')
    
    class Meta:
        abstract = True
    
    def __unicode__(self):
        return '{0} {1}: {2}: {3}-{4}'.format(self.__class__.__name__, self.id,
                                  self.chromosome.name.strip(), 
                                  self.transcription_start, 
                                  self.transcription_end)

        
    @classmethod 
    def add_from_tags(cls,  tag_table):
        connection.close()
        sequencing_run = SequencingRun.objects.get(source_table=tag_table)
        if not sequencing_run.standard: 
            raise Exception('This is not a table marked as "standard," and will not be added to the transcript set.')
        if sequencing_run.type.strip() == 'Gro-Seq':
            cls.add_transcripts_from_groseq(tag_table, sequencing_run)
            
    ################################################
    # Maintenance
    ################################################
    @classmethod
    def force_vacuum(cls):
        '''
        VACUUM ANALYZE all tables.
        '''
        print 'Vacuum analyzing all tables.'
        for model in cls.cell_base.get_transcript_models():
            execute_query_without_transaction('VACUUM ANALYZE "%s";' % (model._meta.db_table))
        
    @classmethod
    def force_vacuum_prep(cls):
        '''
        VACUUM ANALYZE prep tables.
        '''
        print 'Vacuum analyzing prep tables.'
        for model in [cls.cell_base.glass_transcript_prep, cls.cell_base.glass_transcript_source_prep]:
            execute_query_without_transaction('VACUUM ANALYZE "%s";' % (model._meta.db_table))
        multiprocess_all_chromosomes(wrap_force_vacuum, cls)
        
    @classmethod        
    def _force_vacuum(cls, chr_list):
        for chr_id in chr_list:
            print 'Vacuum analyzing chromosome partitions for chr %d' % chr_id
            execute_query_without_transaction('VACUUM ANALYZE "%s_%d";' % (
                                                                    cls.cell_base.glass_transcript_prep._meta.db_table, 
                                                                    chr_id))

            
            
class TranscriptBase(TranscriptionRegionBase):
    '''
    Unique transcribed regions in the genome, as determined via GRO-Seq.
    '''   
    refseq              = models.BooleanField(help_text='Does this transcript overlap with RefSeq transcripts (strand-specific)?')
    
    class Meta:
        abstract = True
        
class GlassTranscriptPrep(TranscriptBase):
    
    class Meta:
        abstract    = True

class GlassTranscript(TranscriptBase):
    distal              = models.BooleanField(help_text='Is this transcript at least 1000 bp away from RefSeq transcripts (not strand-specific)?')
    spliced             = models.NullBooleanField(default=None, help_text='Do we have RNA-Seq confirmation?')
    standard_error      = models.FloatField(null=True, default=None)
    score               = models.FloatField(null=True, default=None)
    
    modified        = models.DateTimeField(auto_now=True)
    created         = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        abstract    = True
          
    ################################################
    # GRO-Seq to transcripts
    ################################################
    @classmethod 
    def add_transcripts_from_groseq(cls,  tag_table, sequencing_run):
        multiprocess_all_chromosomes(wrap_add_transcripts_from_groseq, cls, 
                                     sequencing_run, use_table=tag_table)
         
    @classmethod
    def _add_transcripts_from_groseq(cls, chr_list, sequencing_run):
        for chr_id in chr_list:
            print 'Adding transcripts for chromosome %d' % chr_id
            query = """
                SELECT glass_atlas_%s_%s_prep.save_transcripts_from_sequencing_run(%d, %d,'%s', %d, %d, %d, %d, %d, NULL, NULL);
                """ % (current_settings.GENOME,
                       current_settings.CELL_TYPE.lower(),
                       sequencing_run.id, chr_id, 
                       sequencing_run.source_table.strip(), 
                       MAX_GAP, TAG_EXTENSION, 
                       MAX_EDGE, EDGE_SCALING_FACTOR, DENSITY_MULTIPLIER)
            execute_query(query)
            

    
    ################################################
    # Transcript cleanup and refinement
    ################################################
    @classmethod
    def stitch_together_transcripts(cls, allow_extended_gaps=True, extension_percent='.2', 
                                    set_density=False, null_only=True):
        multiprocess_all_chromosomes(wrap_stitch_together_transcripts, cls, 
                                     allow_extended_gaps, extension_percent, set_density, null_only)
    
    @classmethod
    def _stitch_together_transcripts(cls, chr_list, allow_extended_gaps=True, extension_percent='.2', 
                                     set_density=False, null_only=True):
        '''
        This is tag-level agnostic, stitching based on gap size alone.
        '''
        for chr_id  in chr_list:
            print 'Stitching together transcripts for chromosome %d' % chr_id
            query = """
                SELECT glass_atlas_%s_%s_prep.save_transcripts_from_existing(%d, %d);
                """ % (current_settings.GENOME, 
                       current_settings.CELL_TYPE.lower(),
                       chr_id, MAX_STITCHING_GAP)
            execute_query(query)
            if set_density:
                print 'Setting average tags for preparatory transcripts for chromosome %d' % chr_id
                query = """
                    SELECT glass_atlas_{0}_{1}_prep.set_density({2},{3},{4},{5},{6},{7},{8});
                    """.format(current_settings.GENOME, 
                               current_settings.CELL_TYPE.lower(),
                               chr_id, MAX_EDGE, EDGE_SCALING_FACTOR, 
                               DENSITY_MULTIPLIER, 
                               allow_extended_gaps and 'true' or 'false',
                               extension_percent,
                               null_only and 'true' or 'false')
                execute_query(query)
    
    @classmethod
    def set_density(cls, allow_extended_gaps=True, extension_percent='.2', null_only=True):
        multiprocess_all_chromosomes(wrap_set_density, cls, allow_extended_gaps, extension_percent, null_only)
    
    @classmethod
    def _set_density(cls, chr_list, allow_extended_gaps=True, extension_percent='.2', null_only=True):
        '''
        Force reset average tags for prep DB.
        '''
        for chr_id  in chr_list:
            print 'Setting average tags for preparatory transcripts for chromosome %d' % chr_id
            query = """
                SELECT glass_atlas_{0}_{1}_prep.set_density({2},{3},{4},{5},{6},{7},{8});
                """.format(current_settings.GENOME, 
                           current_settings.CELL_TYPE.lower(),
                           chr_id, MAX_EDGE, EDGE_SCALING_FACTOR, 
                           DENSITY_MULTIPLIER, 
                           allow_extended_gaps and 'true' or 'false',
                           extension_percent,
                           null_only and 'true' or 'false')
            execute_query(query)
            
    @classmethod
    def draw_transcript_edges(cls):
        multiprocess_all_chromosomes(wrap_draw_transcript_edges, cls)
        
    @classmethod
    def _draw_transcript_edges(cls, chr_list):
        for chr_id in chr_list:
            print 'Drawing edges for transcripts for chromosome %d' % chr_id
            query = """
                SELECT glass_atlas_%s_%s%s.draw_transcript_edges(%d);
                """ % (current_settings.GENOME, 
                       current_settings.CELL_TYPE.lower(),
                       current_settings.STAGING,
                       chr_id)
            execute_query(query)           
            
    @classmethod
    def set_scores(cls):
        multiprocess_all_chromosomes(wrap_set_scores, cls)
    
    @classmethod
    def _set_scores(cls, chr_list):
        for chr_id in chr_list:
            print 'Scoring transcripts for chromosome %d' % chr_id
            query = """
                SELECT glass_atlas_%s_%s%s.calculate_scores(%d);
                SELECT glass_atlas_%s_%s%s.calculate_standard_error(%d);
                """ % (current_settings.GENOME,
                       current_settings.CELL_TYPE.lower(),
                       current_settings.STAGING, chr_id,
                       current_settings.GENOME,
                       current_settings.CELL_TYPE.lower(),
                       current_settings.STAGING, chr_id)
            execute_query(query) 

    @classmethod
    def set_start_end_tss(cls):
        multiprocess_all_chromosomes(wrap_set_start_end_tss, cls)
    
    @classmethod
    def _set_start_end_tss(cls, chr_list):
        '''
        We want to use only -1000,+1000 as the TSS for genes,
        but the whole transcript otherwise.
        '''
        schema_name = 'glass_atlas_{0}_{1}'.format(current_settings.GENOME, current_settings.CELL_TYPE.lower())
        for chr_id in chr_list:
            print 'Setting start_end_tss for transcripts for chromosome %d' % chr_id
            
            query = """
            UPDATE {schema_name}.glass_transcript_{chr_id}
                SET start_end_tss = start_end;
                
            UPDATE {schema_name}.glass_transcript_{chr_id} t
                SET start_end_tss = 
                
                CASE WHEN (t.strand = 0) THEN 
                        public.make_box(t.transcription_start - 1000, t.transcription_start + 1000)
                    WHEN (t.strand = 1) THEN 
                        public.make_box(t.transcription_end - 1000, t.transcription_end + 1000)
                    ELSE NULL END
            FROM {schema_name}.glass_transcript_sequence s
            WHERE t.id = s.glass_transcript_id
            AND s.major = true;
            """.format(schema_name=schema_name, chr_id=chr_id)
            execute_query(query)
    
    @classmethod
    def associate_interactions(cls, source_table):
        connection.close()
        sequencing_run = SequencingRun.objects.get(source_table=source_table)
        multiprocess_all_chromosomes(wrap_associate_interactions, cls, sequencing_run)
    
    @classmethod
    def _associate_interactions(cls, chr_list, sequencing_run):
        schema_name = 'glass_atlas_{0}_{1}'.format(current_settings.GENOME, current_settings.CELL_TYPE.lower())
        
        for chr_id in chr_list:
            print 'Associating interactions for chromosome %d' % chr_id
            for strand in (0,1):
                query = """
                    
                    CREATE TABLE {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand}
                        (
                            "chromosome_id" int4 DEFAULT NULL,
                            "glass_transcript_id" int4 DEFAULT NULL,
                            "glass_transcript_2_id" int4 DEFAULT NULL,
                            "sequencing_run_id" int4 DEFAULT NULL,
                            "count" int4 DEFAULT NULL,
                            "chromosome_2_id" int4 DEFAULT NULL,
                            "start_end" box DEFAULT NULL,
                            "strand" int2 DEFAULT NULL
                        );
                    INSERT INTO {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand}
                        (chromosome_id, glass_transcript_id, 
                        sequencing_run_id, "count", chromosome_2_id, start_end, strand) 
                    select {chr_id}, t.id, {sequencing_run_id}, i."count", 
                        i.chromosome_2_id, i.start_end_2, i.strand_2
                    FROM {schema_name}.glass_transcript_{chr_id} t
                    JOIN "{source_table}_{chr_id}_{strand}" i
                    ON t.chromosome_id = i.chromosome_1_id
                    AND t.strand = i.strand_1
                    AND t.start_end_tss && i.start_end_1
                    WHERE t.score >= {min_score}
                    AND t.strand = {strand};
                
                    CREATE INDEX interaction_{chr_id}_{strand}_chr_idx 
                        ON {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand} 
                        USING btree (chromosome_2_id);
                    CREATE INDEX interaction_{chr_id}_{strand}_start_end_idx 
                        ON {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand} 
                        USING gist (start_end);
                    CREATE INDEX interaction_{chr_id}_{strand}_transcript_idx 
                        ON {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand} 
                        USING btree (glass_transcript_id, glass_transcript_2_id);
                    ANALYZE {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand};
                    
                    UPDATE {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand} prep
                    SET glass_transcript_2_id = t2.id
                    FROM {schema_name}.glass_transcript t2
                    WHERE t2.chromosome_id = prep.chromosome_2_id
                    AND t2.strand = prep.strand
                    AND t2.start_end_tss && prep.start_end
                    AND t2.score >= {min_score};
                
                    DELETE FROM {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand} 
                    WHERE glass_transcript_2_id IS NULL;
                    DELETE FROM {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand} 
                    WHERE glass_transcript_id = glass_transcript_2_id;
                    
                    INSERT INTO {schema_name}.glass_transcript_interaction_{chr_id}
                        ("chromosome_id", "glass_transcript_id", 
                        "glass_transcript_2_id", "sequencing_run_id", "count")
                    SELECT * FROM 
                        (SELECT "chromosome_id", "glass_transcript_id", 
                            "glass_transcript_2_id", "sequencing_run_id", SUM("count")
                        FROM {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand}
                        GROUP BY "chromosome_id", "glass_transcript_id", 
                            "glass_transcript_2_id", "sequencing_run_id") der;
                    
                    DROP TABLE {schema_name}.prep_glass_transcript_interaction_{chr_id}_{strand}; 
                    """.format(schema_name=schema_name,
                               source_table=sequencing_run.source_table.strip(),
                               sequencing_run_id=sequencing_run.id,
                               min_score=MIN_SCORE/4,
                               chr_id=chr_id, strand=strand)
                execute_query(query) 
        discard_temp_tables()

    
    @classmethod
    def nearest_genes(cls):
        schema_name = 'glass_atlas_{0}_{1}'.format(current_settings.GENOME, current_settings.CELL_TYPE.lower())
        query = """
        INSERT INTO {schema_name}.nearest_refseq
        (glass_transcript_id, nearest_refseq_transcript_id, distance, minimal_distance) 
        select glass_transcript_id, nearest_refseq_transcript_id, distance, minimal_distance 
        FROM 
         (select 
            DISTINCT nonrefseq.id as glass_transcript_id, 
            refseq.id as nearest_refseq_transcript_id,
            row_number() OVER (PARTITION by nonrefseq.id ORDER BY 
            abs(((refseq.strand = 0)::int*refseq.transcription_start 
                    + (refseq.strand = 1)::int*refseq.transcription_end) -
                ((nonrefseq.strand = 0)::int*nonrefseq.transcription_start 
                        + (nonrefseq.strand = 1)::int*nonrefseq.transcription_end))
            ),
            abs(((refseq.strand = 0)::int*refseq.transcription_start 
                    + (refseq.strand = 1)::int*refseq.transcription_end) -
                ((nonrefseq.strand = 0)::int*nonrefseq.transcription_start 
                        + (nonrefseq.strand = 1)::int*nonrefseq.transcription_end)) as distance,
            least(abs(nonrefseq.transcription_start - refseq.transcription_start), 
                    abs(nonrefseq.transcription_end - refseq.transcription_end),
                    abs(nonrefseq.transcription_start - refseq.transcription_end), 
                    abs(nonrefseq.transcription_end - refseq.transcription_start)) as minimal_distance
            
        FROM {schema_name}.glass_transcript nonrefseq
        
        join {schema_name}.glass_transcript refseq
        on nonrefseq.chromosome_id = refseq.chromosome_id
        
        -- HAS refseq
        join {schema_name}.glass_transcript_sequence seq
        on refseq.id = seq.glass_transcript_id
        
        
        where refseq.score >= {refseq_score}
        and nonrefseq.score >= {nonrefseq_score} 
        and refseq.refseq = true
        and seq.major = true
        and nonrefseq.refseq = false) der
        where der.row_number = 1
        ;  
        """.format(schema_name=schema_name, refseq_score=MIN_SCORE, nonrefseq_score=int(MIN_SCORE/3))
        execute_query(query) 

class FilteredGlassTranscriptManager(models.Manager):
    def get_query_set(self):
        return super(FilteredGlassTranscriptManager, self).get_query_set().filter(score__gte=MIN_SCORE)
    
class FilteredGlassTranscript(object):
    '''
    Acts as a grouping of shared methods for filtered glass transcripts.
    Not a model because the proxy models that will inherit these methods
    complain if they inherit from an abstract method.
    '''
    objects = FilteredGlassTranscriptManager()
    
    ################################################
    # Visualization
    ################################################
    @classmethod
    def generate_bed_file(cls, output_dir):
        '''
        Generates a BED file (http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
        of all the transcripts and their exons for viewing in the UCSC Genome Browser.
        '''
        max_score = cls.objects.aggregate(max=Max('score'))['max']
        transcripts = cls.objects.order_by('chromosome__id','transcription_start')
        
        cls._generate_bed_file(output_dir, transcripts, max_score, strand=0)
        cls._generate_bed_file(output_dir,transcripts, max_score, strand=1)
        
    @classmethod
    def _generate_bed_file(cls, output_dir, transcripts, max_score, strand=0):
        file_name = 'Glass_Transcripts_%s_%d.bed' % (datetime.now().strftime('%Y_%m_%d'), strand)
        file_path = os.path.join(output_dir, file_name)
        
        transcripts = transcripts.filter(strand=strand)
        
        strand_char = strand and '-' or '+'
        output = ('track name=glass_transcripts_{0} description=' \
                    + '"Glass Atlas Transcripts {1} strand" useScore=1 itemRgb=On\n').format(strand, strand_char)
        
        for trans in transcripts:
            # chrom start end name score strand thick_start thick_end colors? exon_count csv_exon_sizes csv_exon_starts
            # Make sure we don't go beyond the length of the chromosome.
            # Transcripts shouldn't due to previous processing, but just in case, we check here as well.
            start = max(0,trans.transcription_start)
            end = min(trans.chromosome.length,trans.transcription_end)
            row = [trans.chromosome.name, str(start), str(end), 
                   'Transcript_' + str(trans.id), str((float(trans.score)/max_score)*1000),
                   trans.strand and '-' or '+', str(start), str(end), 
                   trans.strand and '0,255,0' or '0,0,255']
            
            '''
            # Add in exons
            # Note: 1 bp start and end exons added because UCSC checks that start and end of 
            # whole transcript match the start and end of the first and last blocks, respectively
            exons = cls.cell_base.glass_transcribed_rna.objects.filter(glass_transcript=trans
                                    ).annotate(tags=Sum('glasstranscribedrnasource__tag_count')
                                    ).filter(tags__gt=1).order_by('transcription_start')
            row += [str(exons.count() + 2),
                    ','.join(['1'] +
                             [str(min(trans.transcription_end - ex.transcription_start,
                                ex.transcription_end - ex.transcription_start)) for ex in exons] + ['1']),
                    ','.join(['0'] +
                             [str(max(1, ex.transcription_start - trans.transcription_start)) for ex in exons] 
                                + [str(end - start - 1)]),
                    ]
            '''
            output += '\t'.join(row) + '\n'
        
        f = open(file_path, 'w')
        f.write(output)
        f.close()
        
          
class TranscriptSourceBase(TranscriptModelBase):
    sequencing_run          = models.ForeignKey(SequencingRun)
    tag_count               = models.IntegerField(max_length=12)
    gaps                    = models.IntegerField(max_length=12)
    
    class Meta:
        abstract = True
        
    def __unicode__(self):
        return '%s: %d tags from %s' % (self.__class__.__name__,
                                          self.tag_count,
                                          self.sequencing_run.source_table.strip())
        
class GlassTranscriptSourcePrep(TranscriptSourceBase):
    class Meta:
        abstract = True

class GlassTranscriptSource(TranscriptSourceBase):
    class Meta:
        abstract = True
        
class GlassTranscriptTranscriptionRegionTable(TranscriptModelBase):
    relationship    = models.CharField(max_length=100, choices=[(x,x) 
                                                    for x in ('contains','is contained by','overlaps with','is equal to')])
    
    table_type      = None
    related_class   = None    
    
    class Meta: 
        abstract = True        
    
    def __unicode__(self):
        return 'GlassTranscript %d %s %s' % (self.glass_transcript.id, 
                                             self.relationship.strip(), 
                                             str(self.foreign_key_field()))    
    
    def foreign_key_field(self): 
        return getattr(self, '%s_transcription_region' % self.table_type)
        
class GlassTranscriptSequence(GlassTranscriptTranscriptionRegionTable):
    '''
    Relationship between GlassTranscript and the sequence record it maps to.
    '''
    sequence_transcription_region   = models.ForeignKey(SequenceTranscriptionRegion)
    
    table_type = 'sequence'
    related_class = SequenceTranscriptionRegion
    
    class Meta: 
        abstract = True
           
class GlassTranscriptNonCoding(GlassTranscriptTranscriptionRegionTable):
    '''
    Relationship between GlassTranscript and the non coding region it maps to.
    '''
    non_coding_transcription_region = models.ForeignKey(NonCodingTranscriptionRegion)
    
    table_type = 'non_coding'
    related_class = NonCodingTranscriptionRegion

    class Meta: 
        abstract = True
        
class GlassTranscriptInfrastructure(GlassTranscriptTranscriptionRegionTable):
    '''
    Relationship between GlassTranscript and the infrastructural ncRNA it maps to.
    '''
    infrastructure_transcription_region  = models.ForeignKey(InfrastructureTranscriptionRegion)
    
    table_type      = 'infrastructure'
    related_class   = InfrastructureTranscriptionRegion
    
    class Meta: 
        abstract = True
        
