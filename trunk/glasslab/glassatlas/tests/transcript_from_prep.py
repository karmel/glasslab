'''
Created on Feb 24, 2011

@author: karmel
'''
from __future__ import division
import unittest
from glasslab.glassatlas.tests.base import GlassTestCase
from django.db import connection
from random import randint
from glasslab.sequencing.datatypes.tag import GlassTag
from glasslab.glassatlas.datatypes.transcript import DENSITY_MULTIPLIER,\
    TAG_EXTENSION, MAX_EDGE

class TranscriptPrepTestCase(GlassTestCase):
    ##################################################
    # Associating transcripts
    ##################################################
    
    def test_contains_sequence(self):
        # Transcript gets associated with Serbp1.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 67216950, 67225957
        GlassTag.objects.create(strand=0,
                                chromosome_id=6,
                                start=start, end=end,
                                start_end=(start, 0, end, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 67226000, 67239300
        GlassTag.objects.create(strand=0,
                                chromosome_id=6,
                                start=start_2, end=end_2,
                                start_end=(start_2, 0, end_2, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        # Four versions of Serbp1
        seqs = self.cell_base.glass_transcript_sequence.objects.all()
        self.assertEquals(seqs.count(), 4)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]

        self.assertEquals(seqs[0].relationship, 'contains')
        self.assertEquals(seqs[0].glass_transcript, trans)
        self.assertTrue(seqs.get(sequence_transcription_region__sequence_identifier__sequence_identifier='NM_025814'))
    
    def test_overlaps_with_sequence(self):
        # Transcript gets associated with Serbp1.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 67216950, 67225957
        GlassTag.objects.create(strand=0,
                                chromosome_id=6,
                                start=start, end=end,
                                start_end=(start, 0, end, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 67226000, 67238000
        GlassTag.objects.create(strand=0,
                                chromosome_id=6,
                                start=start_2, end=end_2,
                                start_end=(start_2, 0, end_2, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        connection.close()
        # Four versions of Serbp1
        seqs = self.cell_base.glass_transcript_sequence.objects.all()
        self.cell_base.glass_transcript.draw_transcript_edges()
        self.assertEquals(seqs.count(), 4)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]

        self.assertEquals(seqs[0].relationship, 'overlaps with')
        self.assertEquals(seqs[0].glass_transcript, trans)
        self.assertTrue(seqs.get(sequence_transcription_region__sequence_identifier__sequence_identifier='NM_025814'))
        
    def test_is_contained_by_sequence(self):
        # Transcript gets associated with Serbp1.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 67218000, 67225957
        GlassTag.objects.create(strand=0,
                                chromosome_id=6,
                                start=start, end=end,
                                start_end=(start, 0, end, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 67226000, 67231000
        GlassTag.objects.create(strand=0,
                                chromosome_id=6,
                                start=start_2, end=end_2,
                                start_end=(start_2, 0, end_2, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        # Four versions of Serbp1
        seqs = self.cell_base.glass_transcript_sequence.objects.all()
        self.assertEquals(seqs.count(), 4)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]

        self.assertEquals(seqs[0].relationship, 'is contained by')
        self.assertEquals(seqs[0].glass_transcript, trans)
        self.assertTrue(seqs.get(sequence_transcription_region__sequence_identifier__sequence_identifier='NM_025814'))
     
    def test_is_equal_to_sequence(self):
        # Transcript gets associated with Serbp1.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 67216972, 67225957
        GlassTag.objects.create(strand=0,
                                chromosome_id=6,
                                start=start, end=end,
                                start_end=(start, 0, end, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 67226000, 67239246
        GlassTag.objects.create(strand=0,
                                chromosome_id=6,
                                start=start_2, end=end_2,
                                start_end=(start_2, 0, end_2, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        # Four versions of Serbp1
        seqs = self.cell_base.glass_transcript_sequence.objects.all()
        self.assertEquals(seqs.count(), 4)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]

        self.assertEquals(seqs[0].relationship, 'is equal to')
        self.assertEquals(seqs[0].glass_transcript, trans)
        self.assertTrue(seqs.get(sequence_transcription_region__sequence_identifier__sequence_identifier='NM_025814'))
    
    def test_contains_ncrna(self):
        # Transcript gets associated with ncRNA FR001521: non-protein coding (noncoding) transcript.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 13047530, 13047730
        GlassTag.objects.create(strand=1,
                                chromosome_id=15,
                                start=start, end=end,
                                start_end=(start, 0, end, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 13047740, 13050300
        GlassTag.objects.create(strand=1,
                                chromosome_id=15,
                                start=start_2, end=end_2,
                                start_end=(start_2, 0, end_2, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        
        ncrna = self.cell_base.glass_transcript_non_coding.objects.all()
        self.assertEquals(ncrna.count(), 1)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]

        self.assertEquals(ncrna[0].relationship, 'contains')
        self.assertEquals(ncrna[0].glass_transcript, trans)
        self.assertTrue(ncrna.get(non_coding_transcription_region__non_coding_rna__description__contains='FR001521'))
    
    def test_not_strand_matched(self):
        # Transcript does not get associated with mmu-mir-155.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 84714385, 84714400
        GlassTag.objects.create(strand=1,
                                chromosome_id=16,
                                start=start, end=end,
                                start_end=(start, 0, end, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 84714400, 84714449
        GlassTag.objects.create(strand=1,
                                chromosome_id=16,
                                start=start_2, end=end_2,
                                start_end=(start_2, 0, end_2, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        # Four versions of Serbp1
        seqs = self.cell_base.glass_transcript_sequence.objects.all()
        self.assertEquals(seqs.count(), 0)
        ncrna = self.cell_base.glass_transcript_non_coding.objects.all()
        self.assertEquals(ncrna.count(), 0)
        
    def test_contains_conserved(self):
        # Transcript gets associated with conserved regions.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 133183579, 133183800
        GlassTag.objects.create(strand=0,
                                chromosome_id=1,
                                start=start, end=end,
                                start_end=(start, 0, end, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 133183600, 133185126
        GlassTag.objects.create(strand=0,
                                chromosome_id=1,
                                start=start_2, end=end_2,
                                start_end=(start_2, 0, end_2, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        
        cons = self.cell_base.glass_transcript_conserved.objects.all()
        self.assertEquals(cons.count(), 8)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]

        self.assertEquals(cons[0].relationship, 'contains')
        self.assertEquals(cons[0].glass_transcript, trans)
    
    def test_is_contained_by_patterned(self):
        # Transcript gets associated simple repeat.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 122008915, 122009000
        GlassTag.objects.create(strand=0,
                                chromosome_id=9,
                                start=start, end=end,
                                start_end=(start, 0, end, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 122009005, 122009010
        GlassTag.objects.create(strand=0,
                                chromosome_id=9,
                                start=start_2, end=end_2,
                                start_end=(start_2, 0, end_2, 0)
                                )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        
        patt = self.cell_base.glass_transcript_patterned.objects.all()
        self.assertEquals(patt.count(), 1)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]

        self.assertEquals(patt[0].relationship, 'is contained by')
        self.assertEquals(patt[0].glass_transcript, trans)
        self.assertTrue(patt.get(patterned_transcription_region__type='Simple_repeat'))

    
    ##################################################
    # Calculating scores
    ##################################################
    def test_score_short_transcript(self): 
        # Transcript of length < 1500 should have score = total tag count.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 1000, 1500
        for _ in xrange(0,10):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start + randint(-5,5), end=end + randint(-5,5),
                                    start_end=(start, 0, end, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 1010, 1510
        for _ in xrange(0,20):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start_2 + randint(-5,5), end=end_2 + randint(-5,5),
                                    start_end=(start_2, 0, end_2, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        self.cell_base.glass_transcript.set_scores()
        connection.close()
        self.assertEquals(self.cell_base.glass_transcript.objects.count(), 1)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]
        self.assertEquals(trans.score, 20)
        
    def test_score_1500_transcript(self): 
        # Transcript of length = 1500 should have score = total tag count.
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 1050, 2500
        for _ in xrange(0,25):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start, end=end,
                                    start_end=(start, 0, end, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 1050, 2500
        for _ in xrange(0,20):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start_2, end=end_2,
                                    start_end=(start_2, 0, end_2, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        self.cell_base.glass_transcript.set_scores()
        connection.close()
        self.assertEquals(self.cell_base.glass_transcript.objects.count(), 1)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]
        self.assertEquals(trans.score, 25)
        
    def test_score_long_transcript(self): 
        # Transcript of length > 1500 should have score = total tag count/(lenght/1.5).
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 100053, 103003
        for _ in xrange(0,20):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start, end=end,
                                    start_end=(start, 0, end, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 100053, 103003
        for _ in xrange(0,30):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start_2, end=end_2,
                                    start_end=(start_2, 0, end_2, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        self.cell_base.glass_transcript.set_scores()
        connection.close()
        self.assertEquals(self.cell_base.glass_transcript.objects.count(), 1)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]
        self.assertEquals(trans.score, 15)
    
    ##################################################
    # Calculating average tags
    ##################################################
    def test_transcript_density(self): 
        # Average tags is tags per bp for one run
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 1050, 5500
        for _ in xrange(0,50):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start, end=end,
                                    start_end=(start, 0, end, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        self.assertEquals(self.cell_base.glass_transcript.objects.count(), 1)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]
        self.assertAlmostEquals(trans.density, DENSITY_MULTIPLIER*50/4500)
        
    def test_transcript_density_2(self): 
        # Average tags should be tags per bp per run
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 1050, 2500
        for _ in xrange(0,50):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start, end=end,
                                    start_end=(start, 0, end, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2, end_2 = 1050, 2500
        for _ in xrange(0,20):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start_2, end=end_2,
                                    start_end=(start_2, 0, end_2, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        self.assertEquals(self.cell_base.glass_transcript.objects.count(), 1)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]
        self.assertAlmostEquals(trans.density, DENSITY_MULTIPLIER*70/2/1500)
    
    def test_transcript_start_end_density(self): 
        # Average tags gets transfered over to start_end_density
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 1050, 2000
        for _ in xrange(0,10):
            GlassTag.objects.create(strand=1,
                                    chromosome_id=22,
                                    start=start, end=end,
                                    start_end=(start, 0, end, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        self.assertEquals(self.cell_base.glass_transcript.objects.count(), 1)
        
        trans = self.cell_base.glass_transcript.objects.all()[:1][0]
        self.assertEquals(trans.density, DENSITY_MULTIPLIER*10/1000)
        self.assertEquals(trans.start_end_density, '(%d,%d),(%d,%d)' % (end, trans.density, 
                                                                       start - TAG_EXTENSION, trans.density))
                                                                       
    
    ##################################################
    # Edge creation
    ##################################################
    def test_edge_same_density(self): 
        # With identical densities, transcripts exactly MAX_EDGE apart should unite
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 1550, 2500
        strand, chr_id = 1, 12
        for _ in xrange(0,10):
            GlassTag.objects.create(strand=strand,
                                    chromosome_id=chr_id,
                                    start=start, end=end,
                                    start_end=(start, 0, end, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2 = end + MAX_EDGE + TAG_EXTENSION
        end_2 = start_2 + (end - start)
        
        for _ in xrange(0,10):
            GlassTag.objects.create(strand=strand,
                                    chromosome_id=chr_id,
                                    start=start_2, end=end_2,
                                    start_end=(start_2, 0, end_2, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        self.assertEquals(self.cell_base.glass_transcript.objects.count(), 1)
        
    def test_no_edge_diff_density(self): 
        # With differing densities, transcripts exactly MAX_EDGE apart should unite
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 1550, 2500
        strand, chr_id = 1, 12
        for _ in xrange(0,11):
            GlassTag.objects.create(strand=strand,
                                    chromosome_id=chr_id,
                                    start=start, end=end,
                                    start_end=(start, 0, end, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2 = end + MAX_EDGE + TAG_EXTENSION
        end_2 = start_2 + (end - start)
        
        for _ in xrange(0,10):
            GlassTag.objects.create(strand=strand,
                                    chromosome_id=chr_id,
                                    start=start_2, end=end_2,
                                    start_end=(start_2, 0, end_2, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        self.assertEquals(self.cell_base.glass_transcript.objects.count(), 2)
    
    def test_no_edge_low_density(self): 
        # With low densities, transcripts exactly MAX_EDGE apart should not unite
        self.create_tag_table(sequencing_run_name='sample_run_1', sequencing_run_type='Gro-Seq')
        start, end = 1050, 2500
        strand, chr_id = 1, 12
        for _ in xrange(0,10):
            GlassTag.objects.create(strand=strand,
                                    chromosome_id=chr_id,
                                    start=start, end=end,
                                    start_end=(start, 0, end, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.create_tag_table(sequencing_run_name='sample_run_2', sequencing_run_type='Gro-Seq')
        start_2 = end + MAX_EDGE + TAG_EXTENSION
        end_2 = start_2 + (end - start)
        
        for _ in xrange(0,10):
            GlassTag.objects.create(strand=strand,
                                    chromosome_id=chr_id,
                                    start=start_2, end=end_2,
                                    start_end=(start_2, 0, end_2, 0)
                                    )
        self.cell_base.glass_transcript.add_from_tags(GlassTag._meta.db_table)
        
        self.cell_base.glass_transcript.stitch_together_transcripts()
        self.cell_base.glass_transcript.draw_transcript_edges()
        connection.close()
        self.assertEquals(self.cell_base.glass_transcript.objects.count(), 2)
        
if __name__=='__main__':
    unittest.main()