'''
Created on Nov 12, 2010

@author: karmel

Convenience script for transcript functions.
'''
genome = 'mm9'
cell_type='thiomac'
def sql(genome, cell_type):
    return """
-- Not run from within the codebase, but kept here in case functions need to be recreated.

CREATE OR REPLACE FUNCTION glass_atlas_%s_%s_staging.draw_transcript_edges(chr_id integer)
RETURNS VOID AS $$
DECLARE
    strand integer;
    last_trans glass_atlas_%s_%s_staging.glass_transcript;
    trans glass_atlas_%s_%s_staging.glass_transcript;
    transcript glass_atlas_%s_%s_staging.glass_transcript;
BEGIN
    FOR strand in 0..1
    LOOP
        last_trans := NULL;
        FOR trans IN 
            SELECT * FROM glass_atlas_%s_%s_staging.get_close_transcripts(chr_id, strand)
        LOOP
            -- Initialize the transcript to be saved if necessary
            IF last_trans IS NULL THEN 
                last_trans := trans;
            ELSE
                
                -- Include this tag if it overlaps or has < max_gap bp gap
                -- Else, this is a new transcript; close off the current and restart.
                IF last_trans.transcription_end >= trans.transcription_start THEN 
                    last_trans.transcription_start := (SELECT LEAST(last_trans.transcription_start, trans.transcription_start));
                    last_trans.transcription_end := (SELECT GREATEST(last_trans.transcription_end, trans.transcription_end));
                ELSE
                    transcript := (SELECT glass_atlas_%s_%s_staging.insert_transcript(last_trans));
                    last_trans := trans;
                END IF;
            END IF;
        END LOOP;
        -- And the last one..
        IF last_trans.transcription_start IS NOT NULL THEN
            transcript := (SELECT glass_atlas_%s_%s_staging.insert_transcript(last_trans));
        END IF;
    END LOOP;
        
    PERFORM glass_atlas_%s_%s_staging.insert_transcript_source_records(chr_id);
    PERFORM glass_atlas_%s_%s_staging.insert_associated_transcript_regions(chr_id);
    RETURN;
END;
$$ LANGUAGE 'plpgsql';

CREATE OR REPLACE FUNCTION glass_atlas_%s_%s_staging.get_close_transcripts(chr_id integer, strand integer)
RETURNS SETOF glass_atlas_%s_%s_staging.glass_transcript AS $$
DECLARE
    above_thresh_table text;
    trans record;
    transcript glass_atlas_%s_%s_staging.glass_transcript;    
BEGIN
        -- From prep tables, get close-enough transcripts
    above_thresh_table = 'above_thresh_table_' || chr_id || '_' || strand || '_' || (1000*RANDOM())::int;
        
    EXECUTE 'CREATE TEMP TABLE ' || above_thresh_table  || ' 
        as SELECT t.id, t.strand, t.transcription_start, t.transcription_end, 
                ''0,0,0''::circle as density_circle, point(0,0) as start_density
                FROM  glass_atlas_%s_%s_prep.glass_transcript_' || chr_id || ' t
                JOIN glass_atlas_%s_%s_prep.glass_transcript_source_' || chr_id || ' s
                    ON t.id = s.glass_transcript_id
                WHERE t.strand = ' || strand || '
                GROUP BY t.id, t.strand, t.transcription_start, t.transcription_end
                -- Omit one-tag-wonders and one-run wonders unless they have at least 5 tags
                HAVING avg(s.tag_count) > 1 AND (count(s.sequencing_run_id) > 1 OR avg(s.tag_count) > 5)
            ';
    EXECUTE 'UPDATE ' || above_thresh_table  || ' temp_t
        SET density_circle = t.density_circle, 
            start_density = t.start_density
        FROM glass_atlas_%s_%s_prep.glass_transcript_' || chr_id || ' t
        WHERE t.id = temp_t.id';
    EXECUTE 'CREATE INDEX ' || above_thresh_table || '_density_circle_idx ON ' || above_thresh_table || ' USING gist(density_circle)';
    EXECUTE 'CREATE INDEX ' || above_thresh_table || '_start_density_idx ON ' || above_thresh_table || ' USING gist(start_density)';
        
    FOR trans IN 
        EXECUTE 'SELECT t1.id, t1.strand, t1.transcription_start, t1.transcription_end,
                    GREATEST(max(t2.transcription_end), t1.transcription_end) as max_end
            FROM ' || above_thresh_table  || ' t1
            LEFT OUTER JOIN ' || above_thresh_table  || ' t2 
            ON t1.density_circle @> t2.start_density 
                AND t1.strand = t2.strand 
                AND t1.transcription_start <= t2.transcription_start -- Omit upstream transcripts in the circle 
            GROUP by t1.id, t1.strand, t1.transcription_start, t1.transcription_end
            ORDER by t1.transcription_start ASC'
    LOOP
        IF trans.id IS NOT NULL THEN
            -- Reset transcript record.
            transcript := NULL;
            transcript.chromosome_id = chr_id;
            transcript.strand = trans.strand;
            transcript.transcription_start = trans.transcription_start;
            transcript.transcription_end = trans.max_end;
            
            RETURN NEXT transcript;
        END IF;
    END LOOP;
    RETURN;
END;
$$ LANGUAGE 'plpgsql';

CREATE OR REPLACE FUNCTION glass_atlas_%s_%s_staging.insert_transcript(rec glass_atlas_%s_%s_staging.glass_transcript)
RETURNS glass_atlas_%s_%s_staging.glass_transcript AS $$
DECLARE
    transcript glass_atlas_%s_%s_staging.glass_transcript;
BEGIN
    -- Update record
    EXECUTE 'INSERT INTO glass_atlas_%s_%s_staging.glass_transcript_' || rec.chromosome_id 
    || ' (chromosome_id, strand, transcription_start, transcription_end, '
        || ' start_end, modified, created) '
    || 'VALUES (' || rec.chromosome_id || ',' || rec.strand || ','
        || rec.transcription_start || ',' || rec.transcription_end || ','
        || ' public.make_box(' || rec.transcription_start || ', 0' 
            || ',' ||  rec.transcription_end || ', 0),'
        || ' NOW(), NOW()) RETURNING *' INTO transcript;
    RETURN transcript;
END;
$$ LANGUAGE 'plpgsql';

CREATE OR REPLACE FUNCTION glass_atlas_%s_%s_staging.insert_transcript_source_records(chr_id integer)
RETURNS VOID AS $$
BEGIN 
    EXECUTE 'INSERT INTO glass_atlas_%s_%s_staging.glass_transcript_source_' || chr_id || '
            (chromosome_id, glass_transcript_id, sequencing_run_id, tag_count, gaps)
            SELECT ' || chr_id || ', der.transcript_id, der.sequencing_run_id, 
                SUM(der.tag_count), COUNT(der.glass_transcript_id) - 1
            FROM (
                SELECT trans.id as transcript_id,* 
                FROM glass_atlas_%s_%s_staging.glass_transcript_' || chr_id || ' trans
                JOIN (SELECT * FROM glass_atlas_%s_%s_prep.glass_transcript_' || chr_id || ' t
                    JOIN glass_atlas_%s_%s_prep.glass_transcript_source_' || chr_id || ' s
                    ON t.id = s.glass_transcript_id) source
                ON source.strand = trans.strand
            AND source.start_end <@ trans.start_end) der -- Is contained by
        GROUP BY der.transcript_id, der.sequencing_run_id';
    RETURN;
END;
$$ LANGUAGE 'plpgsql';

CREATE OR REPLACE FUNCTION glass_atlas_%s_%s_staging.insert_associated_transcript_regions(chr_id integer)
RETURNS VOID AS $$
DECLARE
    region_types text[] := ARRAY['sequence','non_coding','conserved','patterned','duped'];
    counter integer;
    table_type text;
BEGIN
    -- Associate any sequencing regions
    FOR counter IN array_lower(region_types,1)..array_upper(region_types,1)
    LOOP
        table_type := region_types[counter];
        EXECUTE 'INSERT INTO glass_atlas_%s_%s_staging.glass_transcript_'
        || table_type || ' (glass_transcript_id, '
        || table_type || '_transcription_region_id, relationship, major)
            (SELECT trans.id, reg.id, 
                (CASE WHEN reg.start_end ~= trans.start_end THEN 
                glass_atlas_%s_%s_staging.glass_transcript_transcription_region_relationship(''is equal to'') 
                WHEN reg.start_end <@ trans.start_end THEN 
                glass_atlas_%s_%s_staging.glass_transcript_transcription_region_relationship(''contains'') 
                WHEN reg.start_end @> trans.start_end THEN 
                glass_atlas_%s_%s_staging.glass_transcript_transcription_region_relationship(''is contained by'') 
                ELSE glass_atlas_%s_%s_staging.glass_transcript_transcription_region_relationship(''overlaps with'') END),
                (CASE WHEN width(reg.start_end # trans.start_end) > width(reg.start_end)::numeric/2 THEN true
                ELSE false END)
            FROM glass_atlas_%s_%s_staging.glass_transcript_' || chr_id || ' trans
            JOIN genome_reference_mm9.' || table_type || '_transcription_region reg
            ON reg.start_end && trans.start_end
            WHERE reg.chromosome_id = ' || chr_id || '
            AND (reg.strand IS NULL OR reg.strand = trans.strand))';
            
    END LOOP;
    
    -- Special case for infrastructure regions
    table_type := 'infrastructure';
    EXECUTE 'INSERT INTO glass_atlas_%s_%s_staging.glass_transcript_'
        || table_type || ' (glass_transcript_id, '
        || table_type || '_transcription_region_id, relationship, major)
            (SELECT trans.id, reg.id, 
                (CASE WHEN reg.start_end ~= trans.start_end THEN 
                glass_atlas_%s_%s_staging.glass_transcript_transcription_region_relationship(''is equal to'') 
                WHEN reg.start_end <@ trans.start_end THEN 
                glass_atlas_%s_%s_staging.glass_transcript_transcription_region_relationship(''contains'') 
                WHEN reg.start_end @> trans.start_end THEN 
                glass_atlas_%s_%s_staging.glass_transcript_transcription_region_relationship(''is contained by'') 
                ELSE glass_atlas_%s_%s_staging.glass_transcript_transcription_region_relationship(''overlaps with'') END),
                (CASE WHEN width(reg.start_end # trans.start_end) > width(reg.start_end)::numeric/2 THEN true
                ELSE false END)
            FROM glass_atlas_%s_%s_staging.glass_transcript_' || chr_id || ' trans
            JOIN genome_reference_mm9.patterned_transcription_region_' || chr_id || ' reg
            ON reg.start_end && trans.start_end
            WHERE reg.type IN (''tRNA'',''rRNA'',''snRNA'',''srpRNA'',''scRNA'',''RNA'') 
            AND (reg.strand IS NULL OR reg.strand = trans.strand))';
    RETURN;
END;
$$ LANGUAGE 'plpgsql';

CREATE OR REPLACE FUNCTION glass_atlas_%s_%s_staging.calculate_scores(chr_id integer)
RETURNS VOID AS $$
DECLARE
    total_runs integer;
BEGIN 
    -- Tag count is scaled avg and max tags: sqrt(avg_tags * max_tags)
    -- Score is tag count divided by the lesser of length/1000 and 2*log(length),
    -- which allows lower tag counts per bp for longer transcripts.
    
    total_runs := (SELECT count(DISTINCT sequencing_run_id) FROM glass_atlas_%s_%s_staging.glass_transcript_source);
	EXECUTE 'UPDATE glass_atlas_%s_%s_staging.glass_transcript_' || chr_id || ' transcript
		    SET score = derived.score 
		    FROM (SELECT 
				    transcript.id,
				    GREATEST(0,
                        SQRT((SUM(source.tag_count)::numeric/' || total_runs || ')
                        *MAX(source.tag_count)::numeric)
                        /LEAST(
                            GREATEST(1000, transcript.transcription_end - transcript.transcription_start + 1)::numeric/1000,
                            2*LOG(transcript.transcription_end - transcript.transcription_start + 1)
                            )
                    ) as score
			    FROM glass_atlas_%s_%s_staging.glass_transcript_' || chr_id || ' transcript 
				JOIN glass_atlas_%s_%s_staging.glass_transcript_source source
			    ON source.glass_transcript_id = transcript.id
			    JOIN glass_atlas_mm9.sequencing_run run
			    ON source.sequencing_run_id = run.id
			    WHERE run.use_for_scoring = true
				    AND transcript.chromosome_id = ' || chr_id || ' 
				    --AND transcript.score IS NULL
			    GROUP BY transcript.id, transcript.transcription_end, transcript.transcription_start
		    ) derived
		    WHERE transcript.id = derived.id';
    EXECUTE 'UPDATE glass_atlas_%s_%s_staging.glass_transcript_' || chr_id || ' transcript
            SET score = 0 WHERE transcript.score IS NULL';
	RETURN;
END;
$$ LANGUAGE 'plpgsql';


CREATE OR REPLACE FUNCTION glass_atlas_%s_%s_staging.calculate_deviation_scores(chr_id integer)
RETURNS VOID AS $$
DECLARE
    total_runs integer;
BEGIN 
    -- We want an easily accessible measure of the expected standard error of this transcript.
    -- Use notx replicates to get a standard deviation of tags, scaled by total tags in each run.
    -- Then we take the sample standard deviation over the square root of the number of samples to get 
    -- standard error of the mean (SEM).
    
    -- This value should be multiplied by total tags in the target run for it to be meaningful!
    
    -- Note that we nest the standard deviation query here three times-- once to ensure that we are using
    -- an array padded out to the desired number of runs (filled with 0),
    -- again to unnest the array so that we can use the Postgresql stddev(expression) function over it,
    -- and again so that we can assign the stddev aggregate to a single update row.
    
    total_runs := (SELECT count(*) FROM glass_atlas_mm9.sequencing_run run
        WHERE run.standard = true
            AND run.type = 'Gro-Seq'
            AND run.wt = true
            AND run.notx = true
            AND run.kla = false
            AND run.other_conditions = false
    );
    EXECUTE 'UPDATE glass_atlas_%s_%s_staging.glass_transcript_' || chr_id || ' transcript
            SET deviation_score = der3.score
            
            FROM (SELECT der2.id, stddev_samp(der2.unnest)/sqrt(' || total_runs || ')::numeric as score 
                
                FROM (SELECT der1.id, unnest(der1.deviation_array) 
                    FROM (SELECT t.id, (array_agg(COALESCE(source_run.tag_count/source_run.total_tags::numeric,0)) 
                            || array_fill(0::numeric,ARRAY[' || total_runs || ']))[1:' || total_runs || '] as deviation_array
                        FROM glass_atlas_%s_%s_staging.glass_transcript_' || chr_id || ' t
                        LEFT OUTER JOIN 
                        (SELECT * FROM glass_atlas_%s_%s_staging.glass_transcript_source_' || chr_id || ' s
                        JOIN glass_atlas_mm9.sequencing_run run
                        ON s.sequencing_run_id = run.id
                        WHERE run.standard = true
                            AND run.type = ''Gro-Seq''
                            AND run.wt = true
                            AND run.notx = true
                            AND run.kla = false
                            AND run.other_conditions = false) source_run
                        ON t.id = source_run.glass_transcript_id
                        GROUP by t.id
                    ) der1
                ) der2
                GROUP BY der2.id
            ) der3
            WHERE transcript.id = der3.id';
            
    RETURN;
END;
$$ LANGUAGE 'plpgsql';

CREATE OR REPLACE FUNCTION glass_atlas_%s_%s_staging.mark_transcripts_spliced(chr_id integer, score_threshold numeric)
RETURNS VOID AS $$
BEGIN
    -- Mark transcripts as spliced if sufficiently high-scoring RNA exists
    EXECUTE 'UPDATE glass_atlas_%s_%s_staging.glass_transcript_' || chr_id || ' transcript 
        SET spliced = true
        FROM glass_atlas_%s_%s_staging.glass_transcribed_rna_' || chr_id || ' transcribed_rna
        WHERE transcribed_rna.glass_transcript_id = transcript.id
            AND transcribed_rna.score >= ' || score_threshold;
    RETURN;
END;
$$ LANGUAGE 'plpgsql'; 


""" % tuple([genome, cell_type]*51)

if __name__ == '__main__':
    print sql(genome, cell_type)
