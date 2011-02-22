'''
Created on Nov 12, 2010

@author: karmel

Convenience script for generated create table statements for transcript tables.
'''

genome = 'no_overlap'
cell_type = 'thiomac'
def sql(genome, cell_type):
    return """
CREATE TABLE "glass_atlas_%s_%s"."glass_transcript" (
    "id" int4 NOT NULL,
    "chromosome_id" int4 DEFAULT NULL,
    "strand" int2 DEFAULT NULL,
    "transcription_start" int8 DEFAULT NULL,
    "transcription_end" int8 DEFAULT NULL,
    "start_end" box DEFAULT NULL,
    "processed" boolean DEFAULT false
);
GRANT ALL ON TABLE "glass_atlas_%s_%s"."glass_transcript" TO  "glass";
CREATE SEQUENCE "glass_atlas_%s_%s"."glass_transcript_id_seq"
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
ALTER SEQUENCE "glass_atlas_%s_%s"."glass_transcript_id_seq" OWNED BY "glass_atlas_%s_%s"."glass_transcript".id;
ALTER TABLE "glass_atlas_%s_%s"."glass_transcript" ALTER COLUMN id SET DEFAULT nextval('"glass_atlas_%s_%s"."glass_transcript_id_seq"'::regclass);
ALTER TABLE ONLY "glass_atlas_%s_%s"."glass_transcript" ADD CONSTRAINT glass_transcript_pkey PRIMARY KEY (id);

-- Tables for each chromosome
CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_1" (
    CHECK ( chromosome_id = 1 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_1_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_1_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (chromosome_id);
CREATE INDEX glass_transcript_1_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (strand);
CREATE INDEX glass_transcript_1_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_2" (
    CHECK ( chromosome_id = 2 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_2_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_2_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_2" USING btree (chromosome_id);
CREATE INDEX glass_transcript_2_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_2" USING btree (strand);
CREATE INDEX glass_transcript_2_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_2" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_3" (
    CHECK ( chromosome_id = 3 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_3_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_3_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_3" USING btree (chromosome_id);
CREATE INDEX glass_transcript_3_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_3" USING btree (strand);
CREATE INDEX glass_transcript_3_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_3" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_4" (
    CHECK ( chromosome_id = 4 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_4_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_4_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_4" USING btree (chromosome_id);
CREATE INDEX glass_transcript_4_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_4" USING btree (strand);
CREATE INDEX glass_transcript_4_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_4" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_5" (
    CHECK ( chromosome_id = 5 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_5_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_5_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_5" USING btree (chromosome_id);
CREATE INDEX glass_transcript_5_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_5" USING btree (strand);
CREATE INDEX glass_transcript_5_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_5" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_6" (
    CHECK ( chromosome_id = 6 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_6_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_6_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_6" USING btree (chromosome_id);
CREATE INDEX glass_transcript_6_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_6" USING btree (strand);
CREATE INDEX glass_transcript_6_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_6" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_7" (
    CHECK ( chromosome_id = 7 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_7_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_7_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_7" USING btree (chromosome_id);
CREATE INDEX glass_transcript_7_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_7" USING btree (strand);
CREATE INDEX glass_transcript_7_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_7" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_8" (
    CHECK ( chromosome_id = 8 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_8_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_8_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_8" USING btree (chromosome_id);
CREATE INDEX glass_transcript_8_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_8" USING btree (strand);
CREATE INDEX glass_transcript_8_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_8" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_9" (
    CHECK ( chromosome_id = 9 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_9_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_9_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_9" USING btree (chromosome_id);
CREATE INDEX glass_transcript_9_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_9" USING btree (strand);
CREATE INDEX glass_transcript_9_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_9" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_10" (
    CHECK ( chromosome_id = 10 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_10_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_10_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_10" USING btree (chromosome_id);
CREATE INDEX glass_transcript_10_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_10" USING btree (strand);
CREATE INDEX glass_transcript_10_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_10" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_11" (
    CHECK ( chromosome_id = 11 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_11_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_11_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_11" USING btree (chromosome_id);
CREATE INDEX glass_transcript_11_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_11" USING btree (strand);
CREATE INDEX glass_transcript_11_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_11" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_12" (
    CHECK ( chromosome_id = 12 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_12_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_12_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_12" USING btree (chromosome_id);
CREATE INDEX glass_transcript_12_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_12" USING btree (strand);
CREATE INDEX glass_transcript_12_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_12" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_13" (
    CHECK ( chromosome_id = 13 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_13_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_13_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_13" USING btree (chromosome_id);
CREATE INDEX glass_transcript_13_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_13" USING btree (strand);
CREATE INDEX glass_transcript_13_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_13" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_14" (
    CHECK ( chromosome_id = 14 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_14_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_14_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_14" USING btree (chromosome_id);
CREATE INDEX glass_transcript_14_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_14" USING btree (strand);
CREATE INDEX glass_transcript_14_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_14" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_15" (
    CHECK ( chromosome_id = 15 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_15_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_15_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_15" USING btree (chromosome_id);
CREATE INDEX glass_transcript_15_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_15" USING btree (strand);
CREATE INDEX glass_transcript_15_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_15" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_16" (
    CHECK ( chromosome_id = 16 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_16_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_16_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_16" USING btree (chromosome_id);
CREATE INDEX glass_transcript_16_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_16" USING btree (strand);
CREATE INDEX glass_transcript_16_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_16" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_17" (
    CHECK ( chromosome_id = 17 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_17_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_17_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_17" USING btree (chromosome_id);
CREATE INDEX glass_transcript_17_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_17" USING btree (strand);
CREATE INDEX glass_transcript_17_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_17" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_18" (
    CHECK ( chromosome_id = 18 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_18_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_18_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_18" USING btree (chromosome_id);
CREATE INDEX glass_transcript_18_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_18" USING btree (strand);
CREATE INDEX glass_transcript_18_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_18" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_19" (
    CHECK ( chromosome_id = 19 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_19_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_19_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_19" USING btree (chromosome_id);
CREATE INDEX glass_transcript_19_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_19" USING btree (strand);
CREATE INDEX glass_transcript_19_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_19" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_20" (
    CHECK ( chromosome_id = 20 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_20_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_20_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_20" USING btree (chromosome_id);
CREATE INDEX glass_transcript_20_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_20" USING btree (strand);
CREATE INDEX glass_transcript_20_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_20" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_21" (
    CHECK ( chromosome_id = 21 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_21_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_21_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_21" USING btree (chromosome_id);
CREATE INDEX glass_transcript_21_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_21" USING btree (strand);
CREATE INDEX glass_transcript_21_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_21" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_22" (
    CHECK ( chromosome_id = 22 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_22_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_22_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_22" USING btree (chromosome_id);
CREATE INDEX glass_transcript_22_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_22" USING btree (strand);
CREATE INDEX glass_transcript_22_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_22" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_23" (
    CHECK ( chromosome_id = 23 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_23_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_23_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_23" USING btree (chromosome_id);
CREATE INDEX glass_transcript_23_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_23" USING btree (strand);
CREATE INDEX glass_transcript_23_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_23" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_24" (
    CHECK ( chromosome_id = 24 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_24_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_24_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_24" USING btree (chromosome_id);
CREATE INDEX glass_transcript_24_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_24" USING btree (strand);
CREATE INDEX glass_transcript_24_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_24" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_25" (
    CHECK ( chromosome_id = 25 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_25_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_25_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_25" USING btree (chromosome_id);
CREATE INDEX glass_transcript_25_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_25" USING btree (strand);
CREATE INDEX glass_transcript_25_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_25" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_26" (
    CHECK ( chromosome_id = 26 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_26_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_26_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_26" USING btree (chromosome_id);
CREATE INDEX glass_transcript_26_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_26" USING btree (strand);
CREATE INDEX glass_transcript_26_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_26" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_27" (
    CHECK ( chromosome_id = 27 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_27_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_27_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_27" USING btree (chromosome_id);
CREATE INDEX glass_transcript_27_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_27" USING btree (strand);
CREATE INDEX glass_transcript_27_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_27" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_28" (
    CHECK ( chromosome_id = 28 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_28_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_28_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_28" USING btree (chromosome_id);
CREATE INDEX glass_transcript_28_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_28" USING btree (strand);
CREATE INDEX glass_transcript_28_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_28" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_29" (
    CHECK ( chromosome_id = 29 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_29_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_29_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_29" USING btree (chromosome_id);
CREATE INDEX glass_transcript_29_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_29" USING btree (strand);
CREATE INDEX glass_transcript_29_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_29" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_30" (
    CHECK ( chromosome_id = 30 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_30_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_30_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_30" USING btree (chromosome_id);
CREATE INDEX glass_transcript_30_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_30" USING btree (strand);
CREATE INDEX glass_transcript_30_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_30" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_31" (
    CHECK ( chromosome_id = 31 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_31_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_31_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_31" USING btree (chromosome_id);
CREATE INDEX glass_transcript_31_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_31" USING btree (strand);
CREATE INDEX glass_transcript_31_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_31" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_32" (
    CHECK ( chromosome_id = 32 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_32_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_32_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_32" USING btree (chromosome_id);
CREATE INDEX glass_transcript_32_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_32" USING btree (strand);
CREATE INDEX glass_transcript_32_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_32" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_33" (
    CHECK ( chromosome_id = 33 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_33_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_33_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_33" USING btree (chromosome_id);
CREATE INDEX glass_transcript_33_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_33" USING btree (strand);
CREATE INDEX glass_transcript_33_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_33" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_34" (
    CHECK ( chromosome_id = 34 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_34_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_34_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_34" USING btree (chromosome_id);
CREATE INDEX glass_transcript_34_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_34" USING btree (strand);
CREATE INDEX glass_transcript_34_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_34" USING gist (start_end);

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_35" (
    CHECK ( chromosome_id = 35 )
) INHERITS ("glass_atlas_%s_%s"."glass_transcript");
CREATE INDEX glass_transcript_35_processed_idx ON "glass_atlas_%s_%s"."glass_transcript_1" USING btree (processed);
CREATE INDEX glass_transcript_35_chr_idx ON "glass_atlas_%s_%s"."glass_transcript_35" USING btree (chromosome_id);
CREATE INDEX glass_transcript_35_strand_idx ON "glass_atlas_%s_%s"."glass_transcript_35" USING btree (strand);
CREATE INDEX glass_transcript_35_start_end_idx ON "glass_atlas_%s_%s"."glass_transcript_35" USING gist (start_end);



CREATE OR REPLACE FUNCTION glass_atlas_%s_%s.glass_transcript_insert_trigger()
RETURNS TRIGGER AS $$
BEGIN
    EXECUTE 'INSERT INTO glass_atlas_%s_%s.glass_transcript_' || NEW.chromosome_id || ' VALUES ('
    || quote_literal(NEW.id) || ','
    || quote_literal(NEW.chromosome_id) || ','
    || quote_literal(NEW.strand) || ','
    || quote_literal(NEW.transcription_start) || ','
    || quote_literal(NEW.transcription_end) || ','
    || 'public.make_box(' || quote_literal(NEW.transcription_start) || ', 0, ' 
        || quote_literal(NEW.transcription_end) || ', 0),'
    || coalesce(quote_literal(NEW.spliced),'NULL') || ','
    || coalesce(quote_literal(NEW.score),'NULL') || ','
    || quote_literal(NEW.modified) || ','
    || quote_literal(NEW.created) || ')';
    RETURN NULL;
END;
$$
LANGUAGE 'plpgsql';

-- Trigger function for inserts on main table
CREATE TRIGGER glass_transcript_trigger
    BEFORE INSERT ON "glass_atlas_%s_%s"."glass_transcript"
    FOR EACH ROW EXECUTE PROCEDURE glass_atlas_%s_%s.glass_transcript_insert_trigger();

CREATE TABLE "glass_atlas_%s_%s"."glass_transcript_source" (
    "id" int4 NOT NULL,
    "glass_transcript_id" int4 DEFAULT NULL,
    "sequencing_run_id" int4 DEFAULT NULL,
    "tag_count" int4 DEFAULT NULL,
    "gaps" int4 DEFAULT NULL
);
GRANT ALL ON TABLE "glass_atlas_%s_%s"."glass_transcript_source" TO  "glass";
CREATE SEQUENCE "glass_atlas_%s_%s"."glass_transcript_source_id_seq"
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
ALTER SEQUENCE "glass_atlas_%s_%s"."glass_transcript_source_id_seq" OWNED BY "glass_atlas_%s_%s"."glass_transcript_source".id;
ALTER TABLE "glass_atlas_%s_%s"."glass_transcript_source" ALTER COLUMN id SET DEFAULT nextval('"glass_atlas_%s_%s"."glass_transcript_source_id_seq"'::regclass);
ALTER TABLE ONLY "glass_atlas_%s_%s"."glass_transcript_source" ADD CONSTRAINT glass_transcript_source_pkey PRIMARY KEY (id);
CREATE INDEX glass_transcript_source_transcript_idx ON "glass_atlas_%s_%s"."glass_transcript_source" USING btree (glass_transcript_id);


""" % tuple([genome, cell_type]*231)

if __name__ == '__main__':
    print sql(genome, cell_type)