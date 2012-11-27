'''
Created on Nov 12, 2012

@author: karmel
'''

from glasslab.config import current_settings
from glasslab.utils.database import SqlGenerator
import csv

class GenomeResourcesSqlGenerator(SqlGenerator):
    ''' Generates the SQL queries for building DB schema. '''
    genome = None

    def __init__(self, genome=None, cell_type=None, staging=None, user=None):
        self.genome = genome or current_settings.GENOME.lower()
        self.staging = staging or current_settings.STAGING
        
        self.schema_name = 'genome_reference_{0}'.format(self.genome)
        super(GenomeResourcesSqlGenerator, self).__init__()

    def all_sql(self):
        s = self.schema()
        s += self.convenience_functions()
        s += self.chromosome()
        s += self.sequence_identifier()
        s += self.sequence_transcription_region()
        s += self.non_coding_rna()
        s += self.non_coding_transcription_region()
        s += self.sequencing_run()
        
        return s
    
    def fill_tables_mm9(self):
        q_set = [self.insert_chromosome_mm9_values(),
                 self.import_ucsc_sequence_mm9_values(),
                 self.insert_sequence_mm9_values(),
                 self.insert_sequence_transcription_region_mm9_values(),
                 self.import_ncrna_org_mm9_values(),
                 self.insert_non_coding_mm9_values(),
                 self.insert_non_coding_transcription_region_mm9_values(),]
         
        return q_set
    
    def cleanup(self):
        return """
        DROP TABLE IF EXISTS "{schema_name}"."refGene";
        DROP TABLE IF EXISTS "{schema_name}"."summary";
        DROP TABLE IF EXISTS "{schema_name}"."mm9_bed";
        """.format(schema_name=self.schema_name)
        
    def schema(self):
        return """
        CREATE SCHEMA "{schema_name}" AUTHORIZATION "{user}";
        """.format(schema_name=self.schema_name, user=self.user)

    def chromosome(self):
        table_name = 'chromosome'
        return """
        CREATE TABLE "{schema_name}"."{table_name}" (
            "id" int4 NOT NULL,
            "name" varchar(25),
            "length" bigint
        );
        """.format(schema_name=self.schema_name, table_name=table_name)\
        + self.pkey_sequence_sql(self.schema_name, table_name)
        
    def sequence_identifier(self):
        table_name = 'sequence_identifier'
        return """
        CREATE TYPE "{schema_name}"."{table_name}_type" AS ENUM('mRNA','RNA');

        CREATE TABLE "{schema_name}"."{table_name}" (
            "id" integer NOT NULL,
            "sequence_identifier" varchar(50),
            "type" "{schema_name}"."{table_name}_type" DEFAULT NULL
        );
        ALTER TABLE ONLY "{schema_name}"."{table_name}" ADD CONSTRAINT {table_name}_sequence_identifier_key UNIQUE (sequence_identifier);
        """.format(schema_name=self.schema_name, table_name=table_name)\
        + self.pkey_sequence_sql(self.schema_name, table_name)
        
    def sequence_transcription_region(self):
        table_name = 'sequence_transcription_region'
        return """
        CREATE TABLE "{schema_name}"."{table_name}" (
            id integer NOT NULL,
            sequence_identifier_id integer,
            chromosome_id integer,
            strand smallint,
            transcription_start bigint,
            transcription_end bigint,
            start_end box,
            start_site_1000 box,
            start_end_1000 box
        );
        CREATE INDEX {table_name}_chr_idx ON "{schema_name}"."{table_name}" USING btree (chromosome_id);
        CREATE INDEX {table_name}_strand_idx ON "{schema_name}"."{table_name}" USING btree (strand);
        CREATE INDEX {table_name}_start_end_idx ON "{schema_name}"."{table_name}" USING gist (start_end);
        CREATE INDEX {table_name}_start_end_1000_idx ON "{schema_name}"."{table_name}" USING gist (start_end_1000);
        """.format(schema_name=self.schema_name, table_name=table_name)\
        + self.pkey_sequence_sql(self.schema_name, table_name)
    
    def non_coding_rna(self):
        table_name = 'non_coding_rna'
        return """
        CREATE TABLE "{schema_name}"."{table_name}" (
            id integer NOT NULL,
            type varchar(100),
            description varchar(255),
            source varchar(255)
        );
        CREATE UNIQUE INDEX {table_name}_type_name_idx ON "{schema_name}"."{table_name}" USING btree ("type","description");
        """.format(schema_name=self.schema_name, table_name=table_name)\
        + self.pkey_sequence_sql(self.schema_name, table_name)
    
    def non_coding_transcription_region(self):
        table_name = 'non_coding_transcription_region'
        return """
        CREATE TABLE "{schema_name}"."{table_name}" (
            id integer NOT NULL,
            non_coding_rna_id integer,
            chromosome_id integer,
            strand smallint,
            transcription_start bigint,
            transcription_end bigint,
            start_end box
        );
        CREATE INDEX {table_name}_chr_idx ON "{schema_name}"."{table_name}" USING btree (chromosome_id);
        CREATE INDEX {table_name}_strand_idx ON "{schema_name}"."{table_name}" USING btree (strand);
        CREATE INDEX {table_name}_start_end_idx ON "{schema_name}"."{table_name}" USING gist (start_end);
        """.format(schema_name=self.schema_name, table_name=table_name)\
        + self.pkey_sequence_sql(self.schema_name, table_name)
         
    def sequencing_run(self):
        table_name = 'sequencing_run'
        return """
        CREATE TABLE "{schema_name}"."{table_name}" (
            "id" int4 NOT NULL,
            "cell_type" varchar(50) DEFAULT NULL,
            "name" varchar(100) DEFAULT NULL,
            "source_table" varchar(100) DEFAULT NULL,
            "total_tags" int8 DEFAULT NULL,
            "modified" timestamp(6) NULL DEFAULT NULL,
            "created" timestamp(6) NULL DEFAULT NULL
        );
        CREATE UNIQUE INDEX {table_name}_source_table_idx ON "{schema_name}"."{table_name}" USING btree (source_table);
        """.format(schema_name=self.schema_name, table_name=table_name) \
        + self.pkey_sequence_sql(schema_name=self.schema_name, table_name=table_name)
    
    def convenience_functions(self):
        return """
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
        
        
    def insert_chromosome_mm9_values(self):
        table_name = 'chromosome'
        return """
        INSERT INTO "{schema_name}"."{table_name}" (name, length) values ('chr1', 197195433), ('chr2', 181748088), 
            ('chr3', 159599784), ('chr4', 155630121), ('chr5', 152537260), ('chr6', 149517038), 
            ('chr7', 152524554), ('chr8', 131738872), ('chr9', 124076173), ('chrX', 166650297), 
            ('chrY', 15902556), ('chr10', 129993256), ('chr11', 121843857), ('chr12', 121257531), 
            ('chr13', 120284313), ('chr14', 125194865), ('chr15', 103494975), ('chr16', 98319151), 
            ('chr17', 95272652), ('chr18', 90772032), ('chr19', 61342431), ('chrM', 16300);
        """.format(schema_name=self.schema_name, table_name=table_name)
    
    def import_ucsc_sequence_mm9_values(self):
        '''
        Create a temp table to be normalized and associated appropriately.
        
        File downloaded at: http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz
        
        This is hardcoded to work with the UCSC download as it is. 
        More flexible import logic can be created here.
        '''
        path_to_file = '../data/refGene.txt'
        f = open(path_to_file)
        output = []
        for l in f:
            fields = l.split('\t')
            output.append("""
                INSERT into "{schema_name}"."refGene" 
                    ("name","chrom","strand","txStart","txEnd") 
                    VALUES ('{0}', '{1}', '{2}', {3}, {4});
                """.format(fields[1], fields[2], fields[3], fields[4], fields[5],
                           schema_name=self.schema_name))
        
        return """
        CREATE TABLE "{schema_name}"."refGene" (
            "name" varchar(50) NOT NULL DEFAULT NULL,
            "chrom" varchar(25) NOT NULL DEFAULT NULL,
            "strand" varchar(1) NOT NULL DEFAULT NULL,
            "txStart" int8 NOT NULL DEFAULT NULL,
            "txEnd" int8 NOT NULL DEFAULT NULL
        );
        """.format(schema_name=self.schema_name) \
        + '\n'.join(output)
    
    def insert_sequence_mm9_values(self):
        table_name = 'sequence_identifier'
        return """
        INSERT INTO "{schema_name}"."{table_name}" ("sequence_identifier", "type") 
            SELECT DISTINCT "name", 
            (CASE WHEN substring("name" from 1 for 2) = 'NM' THEN "{schema_name}".sequence_identifier_type('mRNA')
            WHEN  substring("name" from 1 for 2) = 'NR' THEN "{schema_name}".sequence_identifier_type('RNA')
            ELSE NULL END)
         from "{schema_name}"."refGene";
        """.format(schema_name=self.schema_name, table_name=table_name)
    
    def insert_sequence_transcription_region_mm9_values(self):
        table_name = 'sequence_transcription_region'
        return """
        INSERT INTO "{schema_name}"."{table_name}" 
            (sequence_identifier_id, chromosome_id, strand, transcription_start, transcription_end)
        SELECT 
            seq.id, chr.id,
            (CASE WHEN ref.strand = '-' THEN 1 ELSE 0 END),
            ref."txStart"::bigint, ref."txEnd"::bigint
        FROM "{schema_name}"."sequence_identifier" seq, "{schema_name}"."refGene" ref, 
            "{schema_name}"."chromosome" chr
        WHERE ref."name" = seq.sequence_identifier
            AND ref."chrom" = chr."name";
        
        UPDATE "{schema_name}"."{table_name}"  
        SET start_end = public.make_box("transcription_start",0,"transcription_end",0),
            start_site_1000 = public.make_box(("transcription_start" - 1000),0,("transcription_start" + 1000),0),
            start_end_1000 = public.make_box(("transcription_start" - 1000),0,("transcription_end" + 1000),0)
         WHERE strand = 0;
        
        UPDATE "{schema_name}"."{table_name}"  
        SET start_end = public.make_box("transcription_start",0,"transcription_end",0),
            start_site_1000 = public.make_box(("transcription_end" - 1000),0,("transcription_end" + 1000),0),
            start_end_1000 = public.make_box(("transcription_start" - 1000),0,("transcription_end" + 1000),0)
        WHERE strand = 1;
        """.format(schema_name=self.schema_name, table_name=table_name)

    
    def import_ncrna_org_mm9_values(self):
        '''
        Create a temp table to be normalized and associated appropriately.
        
        File downloaded at: http://www.ncrna.org/frnadb/files/summary.zip
        and: http://www.ncrna.org/frnadb/files/mm9_bed.zip --> mm9_bed/mm9.bed
        
        This is hardcoded to work with the fRNAdb download as it is. 
        More flexible import logic can be created here.
        '''
        path_to_summary = '../data/summary.csv'
        f_summary = csv.reader(open(path_to_summary, 'rb'))
        output = []
        for fields in f_summary:
            # ID,acc,Description,SO name,Oranism,Xref,Length
            if fields[0] == 'ID': continue
            
            if not fields[0] or 'Mus musculus' not in fields[4]: continue
            
            # Shorten and clean up some apostrophes
            fields = [val[:255].replace("'", "''") for val in fields]
            output.append("""
                INSERT into "{schema_name}"."summary" 
                    ("ID","Description","SO name","Xref") 
                    VALUES ('{0}', '{1}', '{2}', '{3}');
                """.format(fields[0], fields[2], fields[3], fields[5],
                           schema_name=self.schema_name))
        
        s = """
        CREATE TABLE "{schema_name}"."summary" (
            "ID" varchar(50) NOT NULL DEFAULT NULL,
            "Description" varchar(255) NOT NULL DEFAULT NULL,
            "SO name" varchar(255) NOT NULL DEFAULT NULL,
            "Xref" varchar(255) NOT NULL DEFAULT NULL
        );
        """.format(schema_name=self.schema_name) \
        + '\n'.join(output)
        
        path_to_bed = '../data/mm9.bed'
        f_bed = open(path_to_bed)
        output = []
        for l in f_bed:
            fields = l.split('\t')
            #chr8    119597933       119597953       FR408228        1000    +       
            output.append("""
                INSERT into "{schema_name}"."mm9_bed" 
                    ("name","chrom","strand","start","end") 
                    VALUES ('{0}', '{1}', '{2}', {3}, {4});
                """.format(fields[3], fields[0], fields[5], fields[1], fields[2],
                           schema_name=self.schema_name))
        
        s += """
        CREATE TABLE "{schema_name}"."mm9_bed" (
            "name" varchar(50) NOT NULL DEFAULT NULL,
            "chrom" varchar(25) NOT NULL DEFAULT NULL,
            "strand" varchar(1) NOT NULL DEFAULT NULL,
            "start" int8 NOT NULL DEFAULT NULL,
            "end" int8 NOT NULL DEFAULT NULL
        );
        """.format(schema_name=self.schema_name) \
        + '\n'.join(output)
        
        print s
        return s
    
    def insert_non_coding_mm9_values(self):
        '''
        fRNAdb from ncrna.org has sequences for all organisms in one large file. 
        We only want those relevant to mm9,
        so we match based on the mm9_bed file provided.
        '''
        table_name = 'non_coding_rna'
        return """
        INSERT INTO "{schema_name}"."{table_name}" ("type", "description", "source") 
        SELECT 
            summary."SO name",
            summary."ID" || ': ' || summary."Description",
            summary."Xref"
        FROM "{schema_name}"."summary" summary
        WHERE
            summary."ID" IN (SELECT "name" from "{schema_name}"."mm9_bed");
        """.format(schema_name=self.schema_name, table_name=table_name)
    
    def insert_non_coding_transcription_region_mm9_values(self):
        '''
        fRNAdb from ncrna.org has sequences for all organisms in one large file. 
        We only want those relevant to mm9,
        so we match based on the mm9_bed file provided.
        '''
        table_name = 'non_coding_transcription_region'
        return """
        INSERT INTO "{schema_name}"."{table_name}" (non_coding_rna_id, chromosome_id, strand, 
        transcription_start, transcription_end, start_end) 
        SELECT 
            nc_rna.id, chr.id, 
            (CASE WHEN mm9_bed.strand = '-' THEN '1' ELSE 0 END),
            mm9_bed.start, mm9_bed."end",
            public.make_box(mm9_bed."start",0,mm9_bed."end",0)
        FROM "{schema_name}"."mm9_bed" mm9_bed, 
            "{schema_name}"."summary" summary, 
            "{schema_name}"."non_coding_rna" nc_rna, 
            "{schema_name}"."chromosome" chr
        WHERE 
            mm9_bed."name" = summary."ID"
            AND nc_rna.description = summary."ID" || ': ' || summary."Description"
            AND chr.name = mm9_bed.chrom;
        """.format(schema_name=self.schema_name, table_name=table_name)