'''
Created on Nov 12, 2010

@author: karmel

Convenience script for generated create table statements for transcribed RNA functions.
'''
genome = 'mm11'
sql = """
-- Not run from within the codebase, but kept here in case functions need to be recreated.
DROP FUNCTION IF EXISTS glass_atlas_%s.update_transcribed_rna_source_records(glass_atlas_%s.glass_transcribed_rna, glass_atlas_%s.glass_transcribed_rna);
DROP FUNCTION IF EXISTS glass_atlas_%s.merge_transcribed_rna(glass_atlas_%s.glass_transcribed_rna, glass_atlas_%s.glass_transcribed_rna);
DROP FUNCTION IF EXISTS glass_atlas_%s.save_transcribed_rna(glass_atlas_%s.glass_transcribed_rna);
DROP FUNCTION IF EXISTS glass_atlas_%s.save_transcribed_rna_from_sequencing_run(integer, integer, text, integer);
DROP FUNCTION IF EXISTS glass_atlas_%s.stitch_transcribed_rna_together(integer, integer);
DROP FUNCTION IF EXISTS glass_atlas_%s.join_overlapping_transcribed_rna(integer);
DROP FUNCTION IF EXISTS glass_atlas_%s.process_transcribed_rna_pair(glass_atlas_%s.glass_transcribed_rna_pair, integer[]);
DROP FUNCTION IF EXISTS glass_atlas_%s.get_overlapping_transcribed_rna(integer);
DROP FUNCTION IF EXISTS glass_atlas_%s.associate_transcribed_rna(integer);
DROP FUNCTION IF EXISTS glass_atlas_%s.mark_transcripts_as_spliced(integer);
DROP TYPE IF EXISTS glass_atlas_%s.glass_transcribed_rna_pair;

CREATE TYPE glass_atlas_%s.glass_transcribed_rna_pair AS ("t1" glass_atlas_%s.glass_transcribed_rna, "t2" glass_atlas_%s.glass_transcribed_rna);

CREATE FUNCTION glass_atlas_%s.update_transcribed_rna_source_records(merged_trans glass_atlas_%s.glass_transcribed_rna, trans glass_atlas_%s.glass_transcribed_rna)
RETURNS VOID AS $$
BEGIN 
	-- Update redundant records: those that already exist for the merge
	UPDATE glass_atlas_%s.glass_transcribed_rna_source merged_assoc
		SET tag_count = (merged_assoc.tag_count + trans_assoc.tag_count),
			gaps = (merged_assoc.gaps + trans_assoc.gaps)
		FROM glass_atlas_%s.glass_transcribed_rna_source trans_assoc
		WHERE merged_assoc.glass_transcribed_rna_id = merged_trans.id
			AND trans_assoc.glass_transcribed_rna_id = trans.id
			AND merged_assoc.sequencing_run_id = trans_assoc.sequencing_run_id;
	-- UPDATE redundant records: those that don't exist for the merge
	UPDATE glass_atlas_%s.glass_transcribed_rna_source
		SET glass_transcribed_rna_id = merged_trans.id
		WHERE glass_transcribed_rna_id = trans.id
		AND sequencing_run_id 
		NOT IN (SELECT sequencing_run_id 
				FROM glass_atlas_%s.glass_transcribed_rna_source
			WHERE glass_transcribed_rna_id = merged_trans.id);
	-- Delete those that remain for the removed transcript.
	DELETE FROM glass_atlas_%s.glass_transcribed_rna_source
		WHERE glass_transcribed_rna_id = trans.id;
	RETURN;
END;
$$ LANGUAGE 'plpgsql';
	
CREATE FUNCTION glass_atlas_%s.merge_transcribed_rna(merged_trans glass_atlas_%s.glass_transcribed_rna, trans glass_atlas_%s.glass_transcribed_rna)
RETURNS glass_atlas_%s.glass_transcribed_rna AS $$
BEGIN
	-- Update the merged transcript
	merged_trans.glass_transcript_id := trans.glass_transcript_id;
	merged_trans.transcription_start := (SELECT LEAST(merged_trans.transcription_start, trans.transcription_start));
	merged_trans.transcription_end := (SELECT GREATEST(merged_trans.transcription_end, trans.transcription_end));
	merged_trans.strand := trans.strand;
	merged_trans.start_end := public.cube(merged_trans.transcription_start, merged_trans.transcription_end);
	RETURN merged_trans;
END;
$$ LANGUAGE 'plpgsql';

CREATE FUNCTION glass_atlas_%s.save_transcribed_rna(rec glass_atlas_%s.glass_transcribed_rna)
RETURNS VOID AS $$
BEGIN 	
	-- Update record
	EXECUTE 'UPDATE glass_atlas_%s.glass_transcribed_rna' 
	|| ' SET'
		|| ' glass_transcript_id = ' || rec.glass_transcript_id || ','
		|| ' strand = ' || rec.strand || ','
		|| ' transcription_start = ' || rec.transcription_start || ','
		|| ' transcription_end = ' || rec.transcription_end || ','
		|| ' start_end = public.cube(' || rec.transcription_start || '::float,' || rec.transcription_end || '::float),'
		|| ' modified = NOW()'
	|| ' WHERE id = ' || rec.id;
	
	RETURN;
END;
$$ LANGUAGE 'plpgsql';

CREATE FUNCTION glass_atlas_%s.save_transcribed_rna_from_sequencing_run(seq_run_id integer, chr_id integer, source_t text, max_gap integer)
RETURNS VOID AS $$
 DECLARE
	rec glass_atlas_%s.glass_transcript_row;
	transcribed_rna glass_atlas_%s.glass_transcribed_rna;
 BEGIN
 	FOR strand IN 0..1 
 	LOOP
		FOR rec IN 
			SELECT * FROM glass_atlas_%s.determine_transcripts_from_sequencing_run(chr_id, strand, source_t, max_gap)
		LOOP		
			-- Saved the Transcribed RNA
			INSERT INTO glass_atlas_%s.glass_transcribed_rna 
				("chromosome_id", "strand", "transcription_start", "transcription_end", "start_end",
				modified, created)
				VALUES (chr_id, strand, rec.transcription_start, rec.transcription_end, 
				public.cube(rec.transcription_start, rec.transcription_end), NOW(), NOW())
				RETURNING * INTO transcribed_rna;
				
			-- Save the record of the sequencing run source
			INSERT INTO glass_atlas_%s.glass_transcribed_rna_source 
				("glass_transcribed_rna_id", "sequencing_run_id", "tag_count", "gaps") 
				VALUES (transcribed_rna.id, seq_run_id, rec.tag_count, rec.gaps);
							
		END LOOP;
	END LOOP;
	
	RETURN;
END;
$$ LANGUAGE 'plpgsql';


CREATE FUNCTION glass_atlas_%s.associate_transcribed_rna(chr_id integer)
RETURNS VOID AS $$
BEGIN
	-- Match transcribed RNA to transcripts, preferring transcripts according to:
	-- 1. Largest transcript containing transcribed RNA
	-- 2. Transcript overlapping with transcribed RNA with minimal distance of extension
	UPDATE glass_atlas_%s.glass_transcribed_rna 
		SET glass_transcript_id = NULL, modified = NOW()
		WHERE chromosome_id = chr_id;
	
	UPDATE glass_atlas_%s.glass_transcribed_rna rna
		SET glass_transcript_id = transcript.id, modified = NOW()
		FROM (SELECT row_number() OVER (
					PARTITION BY rna.id
					ORDER BY (t.transcription_end - t.transcription_start) DESC
				) as row_num,
				rna.id as rna_id, t.id as id
			FROM glass_atlas_%s.glass_transcript t
			JOIN glass_atlas_%s.glass_transcribed_rna rna
			ON t.chromosome_id = rna.chromosome_id
			AND t.strand = rna.strand
			AND t.start_end OPERATOR(public.@>) rna.start_end
			WHERE rna.chromosome_id = chr_id
		) transcript
		WHERE transcript.rna_id = rna.id
		AND transcript.row_num = 1;
	
	UPDATE glass_atlas_%s.glass_transcribed_rna rna
		SET glass_transcript_id = transcript.id, modified = NOW()
		FROM (SELECT row_number() OVER (
					PARTITION BY rna.id
					ORDER BY (GREATEST((t.transcription_start - rna.transcription_start),0)
					+ GREATEST((rna.transcription_end - t.transcription_end),0)) ASC
				) as row_num,
				rna.id as rna_id, t.id as id
			FROM glass_atlas_%s.glass_transcript t
			JOIN glass_atlas_%s.glass_transcribed_rna rna
			ON t.chromosome_id = rna.chromosome_id
			AND t.strand = rna.strand
			AND t.start_end OPERATOR(public.&&) rna.start_end
			WHERE rna.chromosome_id = chr_id
		) transcript
		WHERE transcript.rna_id = rna.id
		AND transcript.row_num = 1;

	RETURN;
END;
$$ LANGUAGE 'plpgsql';

CREATE FUNCTION glass_atlas_%s.stitch_transcribed_rna_together(chr_id integer, allowed_gap integer)
RETURNS VOID AS $$
 DECLARE
	transcribed_rna_group record;
	trans glass_atlas_%s.glass_transcribed_rna;
	merged_trans glass_atlas_%s.glass_transcribed_rna;
	should_merge boolean := false;
	merged boolean := false;
	max_gap integer;
	consumed integer[];
BEGIN
	<<group_loop>>
	FOR transcribed_rna_group IN 
		SELECT
			array_agg(rna.*) as transcribed,
			array_agg(rna.id) as transcribed_ids,
			rna.glass_transcript_id,
			exon.id as exon_id,
			exon.length as exon_length
		FROM "glass_atlas_%s"."glass_transcribed_rna" rna
		LEFT OUTER JOIN (SELECT 
				glass_transcript_id, 
				exon.sequence_transcription_region_id,
				exon.start_end,
				(exon.exon_end - exon.exon_start) as length,
				(reg.transcription_end - reg.transcription_start) as sequence_length,
				exon.id
			FROM "glass_atlas_%s"."glass_transcript_sequence" seq
			JOIN genome_reference_mm9.sequence_exon exon
			ON seq.sequence_transcription_region_id = exon.sequence_transcription_region_id
			JOIN genome_reference_mm9.sequence_transcription_region reg
			ON exon.sequence_transcription_region_id = reg.id
		) exon
		ON rna.glass_transcript_id = exon.glass_transcript_id
		AND rna.start_end operator(public.&&) exon.start_end

		WHERE rna.chromosome_id = chr_id
		GROUP BY rna.glass_transcript_id, exon.id, exon.length, exon.sequence_transcription_region_id, rna.strand, 
			exon.sequence_length
		HAVING count(rna.id) > 1
		ORDER BY exon.sequence_length DESC
	LOOP
		-- Transcribed_rna_group is an array of transcribed RNA and ids of associated glass transcripts and exons
		-- such that the group of transcribed RNA all share the same associated exon id,
		-- and therefore are good candidates for stitching together.
		
		-- BAIL if any of the transcripts have already been consumed.
		IF consumed && transcribed_rna_group.transcribed_ids THEN CONTINUE group_loop;
		END IF;
		
		-- Try to group each transcript, bailing if the gap is passed.
		<<pair_loop>>
		FOR trans IN
			SELECT * FROM unnest(transcribed_rna_group.transcribed)
			ORDER BY start_end ASC
		LOOP
			IF transcribed_rna_group.exon_length IS NULL THEN max_gap := allowed_gap;
			ELSE
				max_gap := (SELECT GREATEST(.5*transcribed_rna_group.exon_length::numeric,allowed_gap))::integer;
			END IF;
			
			IF merged_trans IS NULL THEN merged_trans := trans;
			ELSE
				should_merge := false;
				-- Does this transcript connect?
				IF (merged_trans.transcription_end >= trans.transcription_start) THEN should_merge := true;
				ELSE
					IF (trans.transcription_start - merged_trans.transcription_end) <= max_gap THEN should_merge := true;
					END IF;
				END IF;
				
				IF should_merge = true
					THEN 
						merged_trans := (SELECT glass_atlas_%s.merge_transcribed_rna(merged_trans, trans));
						-- Delete/update old associations
						PERFORM glass_atlas_%s.update_transcribed_rna_source_records(merged_trans, trans);
						EXECUTE 'DELETE FROM glass_atlas_%s.glass_transcribed_rna WHERE id = ' || trans.id;
						merged := true;
						IF (consumed @> ARRAY[merged_trans.id]) = false THEN
							consumed := consumed || merged_trans.id;
						END IF;
						consumed := consumed || trans.id;
				ELSE
					-- We have reached a gap; close off any open merged_transcripts
					IF merged THEN
						PERFORM glass_atlas_%s.save_transcribed_rna(merged_trans);
					END IF;
					-- And reset the merged transcript
					merged_trans := trans;
					merged := false;
				END IF;
			END IF;
		END LOOP;
		
		-- We may have one final merged transcript awaiting a save:
		IF merged THEN
			PERFORM glass_atlas_%s.save_transcribed_rna(merged_trans);
		END IF;
		-- And reset the merged transcript
		merged_trans := NULL;
		merged := false;
	END LOOP;

	RETURN;
END;
$$ LANGUAGE 'plpgsql';

CREATE FUNCTION glass_atlas_%s.join_overlapping_transcribed_rna(chr_id integer)
RETURNS VOID AS $$
 DECLARE
	transcribed_rna_pair glass_atlas_%s.glass_transcribed_rna_pair;
	consumed1 integer[];
BEGIN
	-- Loop until all pairs are consumed
	WHILE (consumed1 IS NULL OR consumed1 > array[]::integer[])
	LOOP
		consumed1 := array[]::integer[];
		FOR transcribed_rna_pair IN
			SELECT * FROM glass_atlas_%s.get_overlapping_transcribed_rna(chr_id)
		LOOP
			-- Transcript 1 overlaps with Transcript 2,
			-- and further has the same glass_transcript_id.
			-- Keep track of consumed transcripts to ensure we do not try to process an already deleted transcript
			consumed1 = (SELECT glass_atlas_%s.process_transcribed_rna_pair(transcribed_rna_pair, consumed1));
		END LOOP;
	END LOOP;

	RETURN;
END;
$$ LANGUAGE 'plpgsql';

CREATE FUNCTION glass_atlas_%s.process_transcribed_rna_pair(transcribed_rna_pair glass_atlas_%s.glass_transcribed_rna_pair, consumed integer[])
RETURNS integer[] AS $$
 DECLARE
	trans glass_atlas_%s.glass_transcribed_rna;
	merged_trans glass_atlas_%s.glass_transcribed_rna;
BEGIN
	merged_trans := transcribed_rna_pair.t1;
	trans := transcribed_rna_pair.t2;
	IF (consumed @> ARRAY[trans.id] = false 
		AND consumed @> ARRAY[merged_trans.id] = false) THEN
			merged_trans := (SELECT glass_atlas_%s.merge_transcribed_rna(merged_trans, trans));
			PERFORM glass_atlas_%s.update_transcribed_rna_source_records(merged_trans, trans);
			EXECUTE 'DELETE FROM glass_atlas_%s.glass_transcribed_rna WHERE id = ' || trans.id;
			PERFORM glass_atlas_%s.save_transcribed_rna(merged_trans);
			consumed = consumed || trans.id;
			consumed = consumed || merged_trans.id;
	END IF;
	RETURN consumed;
END;
$$ LANGUAGE 'plpgsql';

CREATE FUNCTION glass_atlas_%s.get_overlapping_transcribed_rna(chr_id integer)
RETURNS SETOF glass_atlas_%s.glass_transcribed_rna_pair AS $$
BEGIN
	RETURN QUERY
		EXECUTE '(SELECT (transcribed_rna1.*)::glass_atlas_%s.glass_transcribed_rna as t1,'
			|| ' (transcribed_rna2.*)::glass_atlas_%s.glass_transcribed_rna as t2'
		|| ' FROM glass_atlas_%s.glass_transcribed_rna transcribed_rna1'
		|| ' JOIN glass_atlas_%s.glass_transcribed_rna transcribed_rna2'
		|| ' ON transcribed_rna1.chromosome_id = transcribed_rna2.chromosome_id'
		|| ' AND transcribed_rna1.glass_transcript_id = transcribed_rna2.glass_transcript_id'
		|| ' AND transcribed_rna1.strand = transcribed_rna2.strand'
		|| ' AND transcribed_rna1.start_end OPERATOR(public.&&) transcribed_rna2.start_end'
		|| ' AND transcribed_rna1.id != transcribed_rna2.id'
		|| ' WHERE transcribed_rna1.chromosome_id = ' || chr_id
			-- Prevent equivalent run arrays from appearing as two separate rows
			|| ' AND transcribed_rna1.id > transcribed_rna2.id'
		|| ' ORDER BY transcribed_rna1.start_end ASC)';
END;
$$ LANGUAGE 'plpgsql';

CREATE FUNCTION glass_atlas_%s.mark_transcripts_as_spliced(chr_id integer)
RETURNS VOID AS $$
BEGIN
	UPDATE glass_atlas_%s.glass_transcript 
		SET spliced = NULL, modified = NOW()
		WHERE chromosome_id = chr_id;
	UPDATE glass_atlas_%s.glass_transcript 
		SET spliced = true, modified = NOW()
		WHERE chromosome_id = chr_id
			AND id IN (SELECT DISTINCT glass_transcript_id FROM glass_atlas_%s.glass_transcribed_rna);
	RETURN;
END;
$$ LANGUAGE 'plpgsql';

""" % tuple([genome]*81)

print sql
