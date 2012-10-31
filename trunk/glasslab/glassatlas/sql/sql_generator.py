'''
Created on Oct 9, 2012

@author: karmel

'''
from glasslab.config import current_settings
from glasslab.glassatlas.sql import transcripts_from_prep_functions,\
    transcripts_from_tags_functions
from glasslab.utils.database import SqlGenerator

class GlassAtlasSqlGenerator(SqlGenerator):
    ''' Generates the SQL queries for building DB schema. '''
    genome = None
    cell_type = None
    staging = None
    
    def __init__(self, genome=None, cell_type=None, staging=None, user=None):
        self.genome = genome.lower() or current_settings.GENOME.lower()
        self.cell_type = cell_type.lower() or current_settings.CELL_TYPE.lower()
        self.staging = staging.lower() or current_settings.STAGING
        
        self.schema_name_prefix = 'glass_atlas_{0}_{1}'.format(self.genome, self.cell_type)
        super(GlassAtlasSqlGenerator, self).__init__()
    
    def all_sql(self):
        return self.prep_set() + self.final_set()
    
    def prep_set(self):
        prep_suffix = '_prep'
        s = self.schema(prep_suffix)
        s += self.prep_table_main_transcript()
        s += self.table_main_source(prep_suffix)
        
        for chr_id in current_settings.GENOME_CHROMOSOMES:
            s += self.prep_table_chrom_transcript(chr_id)
            s += self.table_chrom_source(chr_id, prep_suffix)
        
        s += self.prep_table_trigger_transcript()
        s += self.table_trigger_source(prep_suffix)
        
        s += self.from_tags_functions()

        return s
            
    def final_set(self):
        s = self.schema(self.staging)
        s += self.table_main_transcript()
        s += self.table_main_source(self.staging)
        for chr_id in current_settings.GENOME_CHROMOSOMES:
            s += self.table_chrom_transcript(chr_id)
            s += self.table_chrom_source(chr_id, self.staging)
        
        s += self.table_trigger_transcript()
        s += self.table_trigger_source(self.staging)
        
        s += self.region_relationship_types()
        s += self.table_region_association('sequence')
        s += self.table_region_association('non_coding')
        s += self.table_region_association('infrastructure')
        
        s += self.peak_features()
        s += self.table_norm_sum()
        
        s += self.from_prep_functions()
        
        return s
    
    # Functions
    def from_tags_functions(self):
        return transcripts_from_tags_functions.sql(self.genome, self.cell_type)
    
    def from_prep_functions(self):
        return transcripts_from_prep_functions.sql(self.genome, self.cell_type, self.staging)
    
    # Tables
    def schema(self, suffix=''):
        return """
        CREATE SCHEMA "{schema_name_prefix}{suffix}" AUTHORIZATION "{user}";
        """.format(schema_name_prefix=self.schema_name_prefix, suffix=suffix, user=self.user)
        
    def prep_table_main_transcript(self):
        table_name = 'glass_transcript'
        return """
        CREATE TABLE "{schema_name_prefix}_prep"."{table_name}" (
            "id" int4 NOT NULL,
            "chromosome_id" int4 DEFAULT NULL,
            "strand" int2 DEFAULT NULL,
            "transcription_start" int8 DEFAULT NULL,
            "transcription_end" int8 DEFAULT NULL,
            "start_end" box DEFAULT NULL,
            "density" float DEFAULT NULL, -- RPKM per run, essentially. But number of bp is configurable.
            "edge" float DEFAULT NULL, -- Length of the allowed edge for joining transcripts.
            "start_density" point DEFAULT NULL,
            "density_circle" circle DEFAULT NULL,
            "refseq" boolean DEFAULT NULL
        );
        """.format(schema_name_prefix=self.schema_name_prefix, table_name=table_name, user=self.user)\
        + self.pkey_sequence_sql(self.schema_name_prefix + '_prep', table_name)
        
    def prep_table_chrom_transcript(self, chr_id):   
        return """ 
        CREATE TABLE "{schema_name_prefix}_prep"."glass_transcript_{chr_id}" (
            CHECK ( chromosome_id = {chr_id} )
        ) INHERITS ("{schema_name_prefix}_prep"."glass_transcript");
        CREATE INDEX glass_transcript_{chr_id}_pkey_idx ON "{schema_name_prefix}_prep"."glass_transcript_{chr_id}" USING btree (id);
        CREATE INDEX glass_transcript_{chr_id}_chr_idx ON "{schema_name_prefix}_prep"."glass_transcript_{chr_id}" USING btree (chromosome_id);
        CREATE INDEX glass_transcript_{chr_id}_strand_idx ON "{schema_name_prefix}_prep"."glass_transcript_{chr_id}" USING btree (strand);
        CREATE INDEX glass_transcript_{chr_id}_start_end_idx ON "{schema_name_prefix}_prep"."glass_transcript_{chr_id}" USING gist (start_end);
        CREATE INDEX glass_transcript_{chr_id}_start_density_idx ON "{schema_name_prefix}_prep"."glass_transcript_{chr_id}" USING gist (start_density);
        CREATE INDEX glass_transcript_{chr_id}_density_circle_idx ON "{schema_name_prefix}_prep"."glass_transcript_{chr_id}" USING gist (density_circle);
        """.format(schema_name_prefix=self.schema_name_prefix, chr_id=chr_id)
        
    def prep_table_trigger_transcript(self):
        return """
        CREATE OR REPLACE FUNCTION {schema_name_prefix}_prep.glass_transcript_insert_trigger()
        RETURNS TRIGGER AS $$
        BEGIN
            EXECUTE 'INSERT INTO {schema_name_prefix}_prep.glass_transcript_' || NEW.chromosome_id || ' VALUES ('
            || quote_literal(NEW.id) || ','
            || quote_literal(NEW.chromosome_id) || ','
            || quote_literal(NEW.strand) || ','
            || quote_literal(NEW.transcription_start) || ','
            || quote_literal(NEW.transcription_end) || ','
            || 'public.make_box(' || quote_literal(NEW.transcription_start) || ', 0, ' 
                || quote_literal(NEW.transcription_end) || ', 0)'
            || '),'
            || 'public.make_box(' || quote_literal((NEW.start_density[1])[0]) || ', ' || quote_literal((NEW.start_density[1])[1]) || ', ' 
                || quote_literal((NEW.start_density[0])[0]) || ', ' || quote_literal((NEW.start_density[0])[1]) || '),'
            || 'circle(' || quote_literal(center(NEW.density_circle)) || ', ' || quote_literal(radius(NEW.density_circle)) || ')'
            || '),'
            || quote_literal(NEW.refseq)
            ;
            RETURN NULL;
        END;
        $$
        LANGUAGE 'plpgsql';
        
        -- Trigger function for inserts on main table
        CREATE TRIGGER glass_transcript_trigger
            BEFORE INSERT ON "{schema_name_prefix}_prep"."glass_transcript"
            FOR EACH ROW EXECUTE PROCEDURE {schema_name_prefix}_prep.glass_transcript_insert_trigger();
        """.format(schema_name_prefix=self.schema_name_prefix)
        
    def table_main_source(self, suffix=''):
        table_name = 'glass_transcript_source'
        return """
        CREATE TABLE "{schema_name_prefix}{suffix}"."" (
            "id" int4 NOT NULL,
            "chromosome_id" int4 DEFAULT NULL,
            "glass_transcript_id" int4 DEFAULT NULL,
            "sequencing_run_id" int4 DEFAULT NULL,
            "tag_count" int4 DEFAULT NULL,
            "gaps" int4 DEFAULT NULL
        );
        """.format(schema_name_prefix=self.schema_name_prefix, table_name=table_name, suffix=suffix)\
        + self.pkey_sequence_sql(self.schema_name_prefix + suffix, table_name)
        
    def table_chrom_source(self, chr_id, suffix=''):
        return """
        CREATE TABLE "{schema_name_prefix}{suffix}"."glass_transcript_source_{chr_id}" (
            CHECK ( chromosome_id = {chr_id} )
        ) INHERITS ("{schema_name_prefix}{suffix}"."glass_transcript_source");
        CREATE INDEX glass_transcript_source_{chr_id}_pkey_idx ON "{schema_name_prefix}{suffix}"."glass_transcript_source_{chr_id}" USING btree (id);
        CREATE INDEX glass_transcript_source_{chr_id}_transcript_idx ON "{schema_name_prefix}{suffix}"."glass_transcript_source_{chr_id}" USING btree (glass_transcript_id);
        CREATE INDEX glass_transcript_source_{chr_id}_sequencing_run_idx ON "{schema_name_prefix}{suffix}"."glass_transcript_source_{chr_id}" USING btree (sequencing_run_id);
        """.format(schema_name_prefix=self.schema_name_prefix, chr_id=chr_id, suffix=suffix)
        
    def table_trigger_source(self, suffix=''):
        return """
        CREATE OR REPLACE FUNCTION {schema_name_prefix}{suffix}.glass_transcript_source_insert_trigger()
        RETURNS TRIGGER AS $$
        BEGIN
            EXECUTE 'INSERT INTO {schema_name_prefix}{suffix}.glass_transcript_source_' || NEW.chromosome_id || ' VALUES ('
            || quote_literal(NEW.id) || ','
            || quote_literal(NEW.chromosome_id) || ','
            || quote_literal(NEW.glass_transcript_id) || ','
            || quote_literal(NEW.sequencing_run_id) || ','
            || quote_literal(NEW.tag_count) || ','
            || quote_literal(NEW.gaps)
            || ')'
            ;
            RETURN NULL;
        END;
        $$
        LANGUAGE 'plpgsql';
        
        -- Trigger function for inserts on main table
        CREATE TRIGGER glass_transcript_source_trigger
            BEFORE INSERT ON "{schema_name_prefix}{suffix}"."glass_transcript_source"
            FOR EACH ROW EXECUTE PROCEDURE {schema_name_prefix}{suffix}.glass_transcript_source_insert_trigger();
        """.format(schema_name_prefix=self.schema_name_prefix, suffix=suffix)
        
    def table_main_transcript(self):
        table_name = 'glass_transcript'
        return """
        CREATE TABLE "{schema_name_prefix}{suffix}"."{table_name}" (
            "id" int4 NOT NULL,
            "chromosome_id" int4 DEFAULT NULL,
            "strand" int2 DEFAULT NULL,
            "transcription_start" int8 DEFAULT NULL,
            "transcription_end" int8 DEFAULT NULL,
            "start_end" box DEFAULT NULL,
            "refseq" boolean DEFAULT NULL,
            "distal" boolean DEFAULT NULL,
            "score" numeric DEFAULT NULL,
            "standard_error" numeric DEFAULT NULL,
            "modified" timestamp(6) NULL DEFAULT NULL,
            "created" timestamp(6) NULL DEFAULT NULL
        );
        """.format(schema_name_prefix=self.schema_name_prefix, table_name=table_name, suffix=self.staging)\
        + self.pkey_sequence_sql(self.schema_name_prefix + self.staging, table_name)
        
    
    def table_chrom_transcript(self, chr_id):
        return """
        CREATE TABLE "{schema_name_prefix}{suffix}"."glass_transcript_{chr_id}" (
            CHECK ( chromosome_id = {chr_id} )
        ) INHERITS ("{schema_name_prefix}{suffix}"."glass_transcript");
        CREATE INDEX glass_transcript_{chr_id}_pkey_idx ON "{schema_name_prefix}{suffix}"."glass_transcript_{chr_id}" USING btree (id);
        CREATE INDEX glass_transcript_{chr_id}_chr_idx ON "{schema_name_prefix}{suffix}"."glass_transcript_{chr_id}" USING btree (chromosome_id);
        CREATE INDEX glass_transcript_{chr_id}_strand_idx ON "{schema_name_prefix}{suffix}"."glass_transcript_{chr_id}" USING btree (strand);
        CREATE INDEX glass_transcript_{chr_id}_score_idx ON "{schema_name_prefix}{suffix}"."glass_transcript_{chr_id}" USING btree (score);
        CREATE INDEX glass_transcript_{chr_id}_distal_idx ON "{schema_name_prefix}{suffix}"."glass_transcript_{chr_id}" USING btree (distal);
        CREATE INDEX glass_transcript_{chr_id}_start_end_idx ON "{schema_name_prefix}{suffix}"."glass_transcript_{chr_id}" USING gist (start_end);
        """.format(schema_name_prefix=self.schema_name_prefix, chr_id=chr_id, suffix=self.staging)
        
    def table_trigger_transcript(self):
        return """
        CREATE OR REPLACE FUNCTION {schema_name_prefix}{suffix}.glass_transcript_insert_trigger()
        RETURNS TRIGGER AS $$
        BEGIN
            EXECUTE 'INSERT INTO {schema_name_prefix}{suffix}.glass_transcript_' || NEW.chromosome_id || ' VALUES ('
            || quote_literal(NEW.id) || ','
            || quote_literal(NEW.chromosome_id) || ','
            || quote_literal(NEW.strand) || ','
            || quote_literal(NEW.transcription_start) || ','
            || quote_literal(NEW.transcription_end) || ','
            || 'public.make_box(' || quote_literal(NEW.transcription_start) || ', 0, ' 
                || quote_literal(NEW.transcription_end) || ', 0),'
            || coalesce(quote_literal(NEW.start_end_density),'NULL') || ','
            || coalesce(quote_literal(NEW.score),'NULL') || ','
            || quote_literal(NEW.modified) || ','
            || quote_literal(NEW.created) || ')';
            RETURN NULL;
        END;
        $$
        LANGUAGE 'plpgsql';
        
        -- Trigger function for inserts on main table
        CREATE TRIGGER glass_transcript_trigger
            BEFORE INSERT ON "{schema_name_prefix}{suffix}"."glass_transcript"
            FOR EACH ROW EXECUTE PROCEDURE {schema_name_prefix}{suffix}.glass_transcript_insert_trigger();

        """.format(schema_name_prefix=self.schema_name_prefix, suffix=self.staging)
        
        
    def region_relationship_types(self):
        return """
        CREATE TYPE "{schema_name_prefix}{suffix}"."glass_transcript_transcription_region_relationship" 
        AS ENUM('contains','is contained by','overlaps with','is equal to');
        """.format(schema_name_prefix=self.schema_name_prefix, suffix=self.staging)
        
    def table_region_association(self, region_type):
        table_name = 'glass_transcript_{type}'.format(type=region_type)
        return """
        CREATE TABLE "{schema_name_prefix}{suffix}"."{table_name}" (
            id integer NOT NULL,
            glass_transcript_id integer,
            {type}_transcription_region_id integer,
            relationship "{schema_name_prefix}{suffix}"."glass_transcript_transcription_region_relationship",
            major boolean
        );
        CREATE INDEX {table_name}_transcript_idx ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree (glass_transcript_id);
        CREATE INDEX {table_name}_major_idx ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree (major);
        """.format(schema_name_prefix=self.schema_name_prefix, table_name=table_name, suffix=self.staging) \
        + self.pkey_sequence_sql(self.schema_name_prefix + self.staging, table_name)
        
        
    def table_norm_sum(self):
        table_name = 'norm_sum'
        return """
        CREATE TABLE "{schema_name_prefix}{suffix}"."{table_name}" (
            "id" int4 NOT NULL,
            "name_1" varchar(100) DEFAULT NULL,
            "name_2" varchar(100) DEFAULT NULL,
            "total_runs_1" int4 DEFAULT NULL,
            "total_runs_2" int4 DEFAULT NULL,
            "total_tags_1" int4 DEFAULT NULL,
            "total_tags_2" int4 DEFAULT NULL,
            "transcript_count" int4 DEFAULT NULL,
            "norm_tags_1" int4 DEFAULT NULL,
            "norm_tags_2" int4 DEFAULT NULL,
            "total_norm_factor" decimal(10,6) DEFAULT NULL,
            "norm_factor" decimal(10,6) DEFAULT NULL,
            "modified" timestamp(6) DEFAULT NULL
        );
        CREATE UNIQUE INDEX "{table_name}_name_idx" ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree(name_1,name_2 ASC NULLS LAST);
        """.format(schema_name_prefix=self.schema_name_prefix, table_name=table_name, suffix=self.staging) \
        + self.pkey_sequence_sql(self.schema_name_prefix + self.staging, table_name)
    
    def nearest_refseq(self):
        table_name = 'nearest_refseq'
        return """
        CREATE TABLE "{schema_name_prefix}{suffix}"."{table_name}" (
            "id" int4 NOT NULL,
            "glass_transcript_id" int4 DEFAULT NULL,
            "nearest_refseq_transcript_id" int4 DEFAULT NULL,
            "distance" int4 DEFAULT NULL, -- Distance from start to start
            "minimal_distance" int4 DEFAULT NULL -- Distance from two closest ends
        );
        CREATE INDEX "{table_name}_transcript_idx" ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree(glass_transcript_id);
        CREATE INDEX "{table_name}_refseq_idx" ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree(nearest_refseq_transcript_id);
        CREATE INDEX "{table_name}_distance_idx" ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree(distance);
        """.format(schema_name_prefix=self.schema_name_prefix, table_name=table_name, suffix=self.staging) \
        + self.pkey_sequence_sql(self.schema_name_prefix + self.staging, table_name)
    
        
    def peak_features(self):
        table_name = 'peak_feature'
        return """
        CREATE TYPE "{schema_name_prefix}{suffix}."glass_transcript_feature_relationship" 
            AS ENUM('contains','is contained by','overlaps with','is equal to','is upstream of','is downstream of');
            
        CREATE TABLE "{schema_name_prefix}{suffix}"."{table_name}" (
            "id" int4 NOT NULL,
            "glass_transcript_id" int4 DEFAULT NULL,
            "glass_peak_id" int4 DEFAULT NULL,
            "sequencing_run_id" int4 DEFAULT NULL,
            "peak_type_id" int4 DEFAULT NULL,
            relationship "glass_atlas_{0}_{1}{suffix}"."glass_transcript_feature_relationship",
            "touches" boolean DEFAULT NULL,
            "length" int4 DEFAULT NULL,
            "tag_count" decimal(8,2) DEFAULT NULL,
            "score" decimal(8,2) DEFAULT NULL,
            "distance_to_tss" int4 DEFAULT NULL
        );
        CREATE INDEX {table_name}_glass_transcript_idx ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree (glass_transcript_id);
        CREATE INDEX {table_name}_peak_type_idx ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree (peak_type_id);
        CREATE INDEX {table_name}_relationship_idx ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree (relationship);
        CREATE INDEX {table_name}_touches_idx ON "{schema_name_prefix}{suffix}"."{table_name}" USING btree (touches);
        """.format(schema_name_prefix=self.schema_name_prefix, table_name=table_name, suffix=self.staging) \
        + self.pkey_sequence_sql(self.schema_name_prefix + self.staging, table_name)
            
    def other_tables(self):
        return """
        CREATE TYPE "glass_atlas_{0}"."sequencing_run_type" AS ENUM('Gro-Seq','RNA-Seq','ChIP-Seq');
        CREATE TABLE "glass_atlas_{0}"."sequencing_run" (
            "id" int4 NOT NULL,
            "type" glass_atlas_{0}.sequencing_run_type DEFAULT NULL,
            "cell_type" varchar(50) DEFAULT NULL,
            "name" varchar(100) DEFAULT NULL,
            "source_table" varchar(100) DEFAULT NULL,
            "description" varchar(255) DEFAULT NULL,
            "total_tags" int8 DEFAULT NULL,
            "percent_mapped" numeric(5,2) DEFAULT NULL,
            "peak_type_id" int4 DEFAULT NULL,
            "standard" boolean DEFAULT false,
            "requires_reload" boolean DEFAULT false,
            "modified" timestamp(6) NULL DEFAULT NULL,
            "created" timestamp(6) NULL DEFAULT NULL
        );
        CREATE UNIQUE INDEX sequencing_run_source_table_idx ON "glass_atlas_{0}"."sequencing_run" USING btree (source_table);
        
        CREATE TABLE "glass_atlas_{0}"."sequencing_run_annotation" (
            "id" int4 NOT NULL,
            "sequencing_run_id" int4 DEFAULT NULL,
            "note" varchar(100) DEFAULT NULL
        );
        CREATE INDEX sequencing_run_annotation_run_idx ON "glass_atlas_{0}"."sequencing_run_annotation" USING btree (sequencing_run_id);
        
        CREATE TABLE "glass_atlas_{0}"."peak_type" (
            "id" int4 NOT NULL,
            "type" varchar(50) DEFAULT NULL,
            "diffuse" boolean DEFAULT NULL
        );
        """.format(self.genome, user=self.user) \
        + self.pkey_sequence_sql('glass_atlas_{0}'.format(self.genome), 'sequencing_run') \
        + self.pkey_sequence_sql('glass_atlas_{0}'.format(self.genome), 'sequencing_run_annotation') \
        + self.pkey_sequence_sql('glass_atlas_{0}'.format(self.genome), 'peak_type')
        
    def convenience_functions(self):
            return """
            CREATE OR REPLACE FUNCTION public.box_equality(left_hand box, right_hand box)
            RETURNS integer AS $$
            BEGIN 
                -- WARNING! Not intended to be an accurate measure of equality!
                -- This exists only so that the Django Admin default DISTINCT queries will work.
                -- For sorting boxes, use gist.
                RETURN -1;
            END;
            $$ LANGUAGE 'plpgsql';
            
            CREATE OPERATOR CLASS public.box_ops
                DEFAULT FOR TYPE box USING btree AS
                    OPERATOR        1       < ,
                    OPERATOR        2       <= ,
                    OPERATOR        3       = ,
                    OPERATOR        4       >= ,
                    OPERATOR        5       > ,
                    FUNCTION        1       public.box_equality(box, box);
            
            CREATE OR REPLACE FUNCTION public.admin_link(id integer)
            RETURNS text AS $$
            BEGIN
                RETURN '<a href="/admin/Transcription_ThioMac/glasstranscriptthiomac/' || id || '" target="_blank">' || id || '</a>';
            END;
            $$ LANGUAGE 'plpgsql';
            
            CREATE OR REPLACE FUNCTION public.admin_link_rna(id integer)
            RETURNS text AS $$
            BEGIN
                RETURN '<a href="/admin/Transcription_ThioMac/glasstranscribedrnathiomac/' || id || '" target="_blank">' || id || '</a>';
            END;
            $$ LANGUAGE 'plpgsql';
            
            
            CREATE OR REPLACE FUNCTION public.ucsc_link_mm9(chr_name character, transcription_start bigint, transcription_end bigint, strand smallint)
            RETURNS text AS $$
            DECLARE
                file_name character(10);
            BEGIN
                IF strand = 0 THEN file_name := 'sense';
                ELSE file_name := 'antisense';
                END IF;
                
                RETURN '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doLoadUrl=submit&amp;hgS_loadUrlName=http%3A%2F%2Fbiowhat.ucsd.edu%2Fkallison%2Fucsc%2Fsessions%2Fglasstranscript_' || file_name ||'_strands.txt&amp;db=mm9&amp;position=' || chr_name || '%3A+' || transcription_start || '-' || transcription_end || '" target="_blank">Strand</a>'
                    || ' | <a href="http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doLoadUrl=submit&amp;hgS_loadUrlName=http%3A%2F%2Fbiowhat.ucsd.edu%2Fkallison%2Fucsc%2Fsessions%2Fall_tracks.txt&amp;db=mm9&amp;position=' || chr_name || '%3A+' || transcription_start || '-' || transcription_end || '" target="_blank">All</a>';
            END;
            $$ LANGUAGE 'plpgsql';
            
            CREATE OR REPLACE FUNCTION public.make_box(x1 numeric, y1 numeric, x2 numeric, y2 numeric)
            RETURNS box AS $$
            DECLARE
                s text;
            BEGIN
                s := '((' || x1 || ',' || y1 || '),(' || x2 || ',' || y2 || '))';
                RETURN s::box;
            END;
            $$ LANGUAGE 'plpgsql';
            
            CREATE OR REPLACE FUNCTION public.make_box(x1 numeric, x2 numeric)
            RETURNS box AS $$
            DECLARE
                s text;
            BEGIN
                s := '((' || x1 || ', 0),(' || x2 || ',0))';
                RETURN s::box;
            END;
            $$ LANGUAGE 'plpgsql';
            """
        
if __name__ == '__main__':
    gen = GlassAtlasSqlGenerator(genome='mm9', cell_type='ThioMac', staging='')
    print gen.nearest_refseq()
