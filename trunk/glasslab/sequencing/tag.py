'''
Created on Sep 27, 2010

@author: karmel
'''
from __future__ import division
from django.db import models, connection
from glasslab.genomereference.datatypes import Chromosome
from glasslab.config import current_settings
from glasslab.utils.datatypes.basic_model import DynamicTable, BoxField
from glasslab.utils.database import execute_query
from glasslab.glassatlas.datatypes.transcript import multiprocess_all_chromosomes,\
    wrap_errors
from glasslab.genomereference.datatypes import SequencingRun

    
def wrap_partition_tables(cls, chr_list): wrap_errors(cls._create_partition_tables, chr_list)
def wrap_translate_from_prep(cls, chr_list): wrap_errors(cls._translate_from_prep, chr_list)
def wrap_set_start_end_box(cls, chr_list): wrap_errors(cls._set_start_end_box, chr_list)
def wrap_set_refseq(cls, chr_list): wrap_errors(cls._set_refseq, chr_list)
def wrap_insert_matching_tags(cls, chr_list): wrap_errors(cls._insert_matching_tags, chr_list)
def wrap_add_indices(cls, chr_list): wrap_errors(cls._add_indices, chr_list)        
    
class GlassTag(DynamicTable):
    '''
    Denormalized version of tag input::
        
        strand    chr    start            sequence_matched
        -       chr5    66529332        ATTAGTTGACAGGAATTAGCTAGGAACCACAGAA

    
    '''
    prep_table = None
    
    chromosome              = models.ForeignKey(Chromosome)
    strand                  = models.IntegerField(max_length=1)
    start                   = models.IntegerField(max_length=12)
    end                     = models.IntegerField(max_length=12)
    
    start_end               = BoxField(max_length=255, help_text='This is a placeholder for the PostgreSQL box type.', null=True)
    
    refseq                  = models.NullBooleanField(default=None)
    
     
    @classmethod        
    def create_prep_table(cls, name):
        '''
        Create table that will be used for to upload preparatory tag information,
        dynamically named.
        '''
        cls.set_prep_table('prep_%s' % name)
        table_sql = """
        CREATE TABLE "%s" (
            strand_char character(1) default NULL,
            chromosome varchar(20),
            "start" bigint,
            sequence_matched varchar(100)
        );
        """ % cls.prep_table
        execute_query(table_sql)
    
    @classmethod 
    def set_prep_table(cls, name): 
        cls.prep_table = '%s"."%s' % (current_settings.CURRENT_SCHEMA, name)
        
    @classmethod                
    def create_parent_table(cls, name):
        '''
        Create table that will be used for these tags,
        dynamically named.
        '''
        name = 'tag_%s' % name
        cls.set_table_name(name)
        table_sql = """
        CREATE TABLE "%s" (
            id integer NOT NULL,
            chromosome_id integer default NULL,
            strand smallint default NULL,
            "start" bigint,
            "end" bigint,
            start_end box,
            refseq boolean default false
        );
        CREATE SEQUENCE "%s_id_seq"
            START WITH 1
            INCREMENT BY 1
            NO MINVALUE
            NO MAXVALUE
            CACHE 1;
        ALTER SEQUENCE "%s_id_seq" OWNED BY "%s".id;
        ALTER TABLE "%s" ALTER COLUMN id SET DEFAULT nextval('"%s_id_seq"'::regclass);
        ALTER TABLE ONLY "%s" ADD CONSTRAINT %s_pkey PRIMARY KEY (id);
        """ % (cls._meta.db_table, 
               cls._meta.db_table,
               cls._meta.db_table, cls._meta.db_table,
               cls._meta.db_table, cls._meta.db_table,
               cls._meta.db_table, cls.name)
        execute_query(table_sql)
        
        cls.table_created = True
    
    @classmethod
    def create_partition_tables(cls):
        '''
        Can't be multiprocessed; too many attempts to ANALYZE at once.
        '''
        
        for chr_id in current_settings.GENOME_CHROMOSOMES:
            table_sql = """
            CREATE TABLE "%s_%d" (
                CHECK ( chromosome_id = %d )
            ) INHERITS ("%s");""" % (cls._meta.db_table,
                                     chr_id, chr_id,
                                     cls._meta.db_table)
            execute_query(table_sql)
        
        trigger_sql = '''
            CREATE OR REPLACE FUNCTION %s.glass_tag_insert_trigger()
            RETURNS TRIGGER AS $$
            DECLARE
                refseq_string text;
                start_end_string text;
            BEGIN
                IF NEW.refseq IS NULL THEN
                    refseq_string := 'NULL';
                ELSE
                    refseq_string := NEW.refseq::text;
                END IF;
                
                IF NEW.start_end IS NULL THEN
                    start_end_string := 'NULL';
                ELSE
                    start_end_string := 'public.make_box(' || quote_literal(NEW.start) || ', 0, ' 
                    || quote_literal(NEW."end") || ', 0)';
                END IF;
                EXECUTE 'INSERT INTO "%s_' || NEW.chromosome_id || '" VALUES ('
                || quote_literal(NEW.id) || ','
                || quote_literal(NEW.chromosome_id) || ','
                || quote_literal(NEW.strand) || ','
                || quote_literal(NEW.start) || ','
                || quote_literal(NEW."end") || ','
                || start_end_string || ','
                || refseq_string || ')';
                RETURN NULL;
            END;
            $$
            LANGUAGE 'plpgsql';
            
            -- Trigger function for inserts on main table
            CREATE TRIGGER glass_tag_trigger
                BEFORE INSERT ON "%s"
                FOR EACH ROW EXECUTE PROCEDURE %s.glass_tag_insert_trigger();
        ''' % (current_settings.CURRENT_SCHEMA,
               cls._meta.db_table,
               cls._meta.db_table,
               current_settings.CURRENT_SCHEMA)
        execute_query(trigger_sql)
        
    @classmethod
    def associate_chromosome(cls):
        '''
        Denormalize chromosome value.
        '''
        update_query = """
        UPDATE "%s" tag SET chromosome_id = chr.id
        FROM "%s" chr WHERE tag.chromosome = chr.name;
        """ % (cls._meta.db_table, Chromosome._meta.db_table,)
        execute_query(update_query)
    
    @classmethod
    def translate_from_prep(cls):
        multiprocess_all_chromosomes(wrap_translate_from_prep, cls)
        
    @classmethod
    def _translate_from_prep(cls, chr_list):
        '''
        Take uploaded values and streamline into ints and sequence ends.
        '''
        for chr_id in chr_list:
            update_query = """
            INSERT INTO "%s_%d" (chromosome_id, strand, "start", "end", start_end, refseq)
            SELECT * FROM (
                SELECT %d, (CASE WHEN prep.strand_char = '-' THEN 1 ELSE 0 END), 
                prep."start", (prep."start" + char_length(prep.sequence_matched)),
                public.make_box(prep."start", 0, (prep."start" + char_length(prep.sequence_matched)), 0),
                NULL::boolean
            FROM "%s" prep
            JOIN "%s" chr ON chr.name = prep.chromosome
            WHERE chr.id = %d) derived;
            """ % (cls._meta.db_table,chr_id,
                   chr_id, cls.prep_table,
                   Chromosome._meta.db_table,
                   chr_id)
            execute_query(update_query)
                        
    @classmethod
    def set_start_end_box(cls):
        multiprocess_all_chromosomes(wrap_set_start_end_box, cls)
        
    @classmethod
    def _set_start_end_box(cls, chr_list):
        '''
        Create type box field for faster interval searching with the PostgreSQL box.
        '''
        for chr_id in chr_list:
            update_query = """
            UPDATE "%s" set start_end = public.make_box("start", 0, "end", 0) WHERE chromosome_id = %d;
            """ % (cls._meta.db_table, chr_id)
            execute_query(update_query)
    
    @classmethod
    def set_refseq(cls):
        multiprocess_all_chromosomes(wrap_set_refseq, cls)
        
    @classmethod
    def _set_refseq(cls, chr_list):
        '''
        Create type box field for faster interval searching with the PostgreSQL box.
        '''
        for chr_id in chr_list:
            print 'Setting Refseq status for chromosome {0}'.format(chr_id)
            update_query = """
            UPDATE "{0}_{1}" tag 
            SET refseq = NULL; 

            UPDATE "{0}_{1}" tag 
            SET refseq = true 
            FROM glass_atlas_{2}_refseq.glass_transcript_{1} ref
            WHERE ref.start_end && tag.start_end
            AND ref.strand = tag.strand
            AND tag.refseq IS NULL;

            UPDATE "{0}_{1}" tag 
            SET refseq = false 
            WHERE refseq IS NULL;
            """.format(cls._meta.db_table, chr_id, current_settings.GENOME)
            execute_query(update_query)
            
    @classmethod
    def delete_prep_columns(cls):
        '''
        .. warning:: DELETES the associated sequence and strand_char columns to conserve space.
        
        '''
        table_sql = """ 
        ALTER TABLE "%s" DROP COLUMN sequence_matched;
        ALTER TABLE "%s" DROP COLUMN strand_char;
        ALTER TABLE "%s" DROP COLUMN chromosome;
        """ % (cls._meta.db_table, cls._meta.db_table, cls._meta.db_table)
        execute_query(table_sql)
        
        cls.table_created = False

    @classmethod
    def add_indices(cls):
        multiprocess_all_chromosomes(wrap_add_indices, cls)
        
    @classmethod
    def _add_indices(cls, chr_list):
        for chr_id in chr_list:
            update_query = """
            CREATE INDEX %s_%d_pkey_idx ON "%s_%d" USING btree (id);
            CREATE INDEX %s_%d_chr_idx ON "%s_%d" USING btree (chromosome_id);
            CREATE INDEX %s_%d_strand_idx ON "%s_%d" USING btree (strand);
            CREATE INDEX %s_%d_start_end_idx ON "%s_%d" USING gist (start_end);
            ANALYZE "%s_%d";
            """ % (cls.name, chr_id, cls._meta.db_table, chr_id,
                   cls.name, chr_id, cls._meta.db_table, chr_id,
                   cls.name, chr_id, cls._meta.db_table, chr_id,
                   cls.name, chr_id, cls._meta.db_table, chr_id,
                   cls._meta.db_table, chr_id)
            execute_query(update_query)

    @classmethod 
    def add_record_of_tags(cls):
        '''
        Add SequencingRun record with the details of this run.
        
        Should be called only after all tags have been added.
        '''
        connection.close()
        # If possible, retrieve bowtie stats
        s, created = SequencingRun.objects.get_or_create(source_table=cls._meta.db_table)
        return s