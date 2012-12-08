'''
Created on Dec 7, 2012

@author: karmel

HiC gives us data about interaction sites in the genome.
We want to capture that information and use it to identify
transcripts that interact.
'''
from __future__ import division
from django.db import models, connection
from glasslab.utils.datatypes.genome_reference import Chromosome
from glasslab.config import current_settings
from glasslab.utils.datatypes.basic_model import BoxField
from glasslab.utils.database import execute_query
from glasslab.glassatlas.datatypes.transcript import multiprocess_all_chromosomes
from glasslab.glassatlas.datatypes.metadata import SequencingRun
import re
from glasslab.sequencing.datatypes.tag import GlassSequencingOutput,\
    wrap_translate_from_prep, wrap_add_indices

class GlassInteraction(GlassSequencingOutput):
    '''
    MakeTagDirectory HiC output::
        
        chr_1    start_1    strand_1    count    length_1    chr_2    start_2    strand_2    length_2
        chrY    2902017    0    1.0    57    chr11    39571330    0    57

    
    '''
    prep_table = None
    
    chromosome_1            = models.ForeignKey(Chromosome)
    strand_1                = models.IntegerField(max_length=1)
    start_1                 = models.IntegerField(max_length=12)
    end_1                   = models.IntegerField(max_length=12)
    
    start_end_1             = BoxField(max_length=255, help_text='This is a placeholder for the PostgreSQL box type.', null=True)
    
    chromosome_2            = models.ForeignKey(Chromosome)
    strand_2                = models.IntegerField(max_length=1)
    start_2                 = models.IntegerField(max_length=12)
    end_2                   = models.IntegerField(max_length=12)
    
    start_end_2             = BoxField(max_length=255, help_text='This is a placeholder for the PostgreSQL box type.', null=True)
    
    count                   = models.IntegerField(max_length=12)
     
    @classmethod        
    def create_prep_table(cls, name):
        '''
        Create table that will be used for to upload preparatory tag information,
        dynamically named.
        '''
        cls.set_prep_table('prep_%s' % name)
        table_sql = """
        CREATE TABLE "%s" (
            chromosome_1 varchar(20),
            "start_1" bigint,
            strand_1 int2,
            count int4,
            length_1 int4,
            chromosome_2 varchar(20),
            "start_2" bigint,
            strand_2 int2,
            length_2 int4
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
        name = 'interaction__%s' % name
        cls.set_table_name(name)
        table_sql = """
        CREATE TABLE "%s" (
            id integer NOT NULL,
            chromosome_1_id integer default NULL,
            strand_1 smallint default NULL,
            "start_1" bigint,
            "end_1" bigint,
            start_end_1 box,
            
            chromosome_2_id integer default NULL,
            strand_2 smallint default NULL,
            "start_2" bigint,
            "end_2" bigint,
            start_end_2 box,
            
            count int4
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
                CHECK ( chromosome_1_id = %d )
            ) INHERITS ("%s");""" % (cls._meta.db_table,
                                     chr_id, chr_id,
                                     cls._meta.db_table)
            execute_query(table_sql)
        
        trigger_sql = '''
            CREATE OR REPLACE FUNCTION %s.glass_interaction_insert_trigger()
            RETURNS TRIGGER AS $$
            DECLARE
                start_end_1_string text;
                start_end_2_string text;
            BEGIN
                IF NEW.start_end_1 IS NULL THEN
                    start_end_1_string := 'NULL';
                ELSE
                    start_end_1_string := 'public.make_box(' || quote_literal(NEW.start_1) || ', 0, ' 
                    || quote_literal(NEW."end_1") || ', 0)';
                END IF;
                IF NEW.start_end_2 IS NULL THEN
                    start_end_2_string := 'NULL';
                ELSE
                    start_end_2_string := 'public.make_box(' || quote_literal(NEW.start_2) || ', 0, ' 
                    || quote_literal(NEW."end_2") || ', 0)';
                END IF;
                EXECUTE 'INSERT INTO "%s_' || NEW.chromosome_1_id || '" VALUES ('
                || quote_literal(NEW.id) || ','
                || quote_literal(NEW.chromosome_1_id) || ','
                || quote_literal(NEW.strand_1) || ','
                || quote_literal(NEW.start_1) || ','
                || quote_literal(NEW."end_1") || ','
                || start_end_1_string || ','
                || quote_literal(NEW.chromosome_2_id) || ','
                || quote_literal(NEW.strand_2) || ','
                || quote_literal(NEW.start_2) || ','
                || quote_literal(NEW."end_2") || ','
                || start_end_2_string || ','
                || quote_literal(NEW."count") || ')';
                RETURN NULL;
            END;
            $$
            LANGUAGE 'plpgsql';
            
            -- Trigger function for inserts on main table
            CREATE TRIGGER glass_interaction_trigger
                BEFORE INSERT ON "%s"
                FOR EACH ROW EXECUTE PROCEDURE %s.glass_interaction_insert_trigger();
        ''' % (current_settings.CURRENT_SCHEMA,
               cls._meta.db_table,
               cls._meta.db_table,
               current_settings.CURRENT_SCHEMA)
        execute_query(trigger_sql)
        
    @classmethod
    def translate_from_prep(cls):
        multiprocess_all_chromosomes(wrap_translate_from_prep, cls)
        
    @classmethod
    def _translate_from_prep(cls, chr_list):
        '''
        Take uploaded values and streamline into ints and sequence ends.
        
        Note that starts are 1-indexed, so are converted here to 0-index.
        Also, length is inclusive, so subtract 1.
        '''
        for chr_id in chr_list:
            update_query = """
            INSERT INTO "{0}_{chr_id}" (chromosome_1_id, strand_1, "start_1", "end_1", start_end_1,
                chromosome_2_id, strand_2, "start_2", "end_2", start_end_2, count)
            SELECT * FROM (
                SELECT {chr_id}, prep.strand_1, 
                prep."start_1"-1, (prep."start_1"-1 + prep.length_1-1),
                public.make_box(prep."start_1"-1, 0, (prep."start_1"-1 + prep.length_1-1), 0),
                SELECT chr_2.id, prep.strand_2, 
                prep."start_2"-1, (prep."start_2"-1 + prep.length_2-1),
                public.make_box(prep."start_2"-1, 0, (prep."start_2"-1 + prep.length_2-1), 0),
                prep.count
            FROM "{1}" prep
            JOIN "{2}" chr_1 ON chr_1.name = prep.chromosome_1
            JOIN "{2}" chr_2 ON chr_2.name = prep.chromosome_2
            WHERE chr_1.id = {chr_id}) derived;
            """.format(cls._meta.db_table, cls.prep_table,
                   Chromosome._meta.db_table,
                   chr_id=chr_id)
            execute_query(update_query)
                        
    @classmethod
    def add_indices(cls):
        multiprocess_all_chromosomes(wrap_add_indices, cls)
        
    @classmethod
    def _add_indices(cls, chr_list):
        for chr_id in chr_list:
            update_query = """
            CREATE INDEX {0}_{chr_id}_pkey_idx ON "{1}_{chr_id}" USING btree (id);
            CREATE INDEX {0}_{chr_id}_chr_1_idx ON "{1}_{chr_id}" USING btree (chromosome_1_id);
            CREATE INDEX {0}_{chr_id}_strand_1_idx ON "{1}_{chr_id}" USING btree (strand_1);
            CREATE INDEX {0}_{chr_id}_start_end_1_idx ON "{1}_{chr_id}" USING gist (start_end_1);
            CREATE INDEX {0}_{chr_id}_chr_2_idx ON "{1}_{chr_id}" USING btree (chromosome_2_id);
            CREATE INDEX {0}_{chr_id}_strand_2_idx ON "{1}_{chr_id}" USING btree (strand_2);
            CREATE INDEX {0}_{chr_id}_start_end_2_idx ON "{1}_{chr_id}" USING gist (start_end_2);
            ANALYZE "{1}_{chr_id}";
            """.format(cls.name, cls._meta.db_table, chr_id=chr_id)
            execute_query(update_query)

    @classmethod 
    def add_record_of_tags(cls, description='', type='HiC', standard=False, stats_file=None):
        '''
        Add SequencingRun record with the details of this run.
        
        Should be called only after all tags have been added.
        '''
        connection.close()
        # If possible, retrieve bowtie stats
        wt, notx, kla, other_conditions, timepoint = cls.parse_attributes_from_name()
        s, created = SequencingRun.objects.get_or_create(source_table=cls._meta.db_table,
                                        defaults={'name': cls.name, 
                                                  'total_tags': cls.objects.count(),
                                                  'description': description,
                                                  'cell_type': current_settings.CELL_TYPE,
                                                  'type': type,
                                                  'standard': standard,
                                                  'wt': wt,
                                                  'notx': notx,
                                                  'kla': kla,
                                                  'other_conditions': other_conditions,
                                                  'timepoint': timepoint,
                                                  'strain': wt and 'C57Bl6' or None,
                                                   }
                                               )
        if not created: 
            s.total_tags = cls.objects.count()
            s.timepoint = timepoint 
            s.save() 
        return s

