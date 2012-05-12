'''
Created on Dec 22, 2010

@author: karmel
'''
from glasslab.glassatlas.datatypes.transcript import GlassTranscript,\
    FilteredGlassTranscript, GlassTranscriptNucleotides, GlassTranscriptSource,\
    GlassTranscriptSequence, GlassTranscriptConserved, GlassTranscriptPatterned,\
    GlassTranscriptNonCoding, CellTypeBase, GlassTranscriptSourcePrep,\
    GlassTranscriptPrep, GlassTranscriptInfrastructure, GlassTranscriptDuped
from glasslab.config import current_settings
from glasslab.glassatlas.datatypes.transcribed_rna import GlassTranscribedRna,\
    GlassTranscribedRnaSource
from django.db import models
from glasslab.glassatlas.datatypes.feature import PeakFeature
from glasslab.glassatlas.datatypes.metadata import TranscriptClass
from glasslab.glassatlas.datatypes.label import GlassTranscriptLabel

CELL_TYPE = 'BMDC'

class BMDCBase(CellTypeBase):
    cell_type = CELL_TYPE
           
    @property
    def glass_transcript(self): return GlassTranscriptBMDC
    @property
    def glass_transcript_prep(self): return GlassTranscriptPrepBMDC
    @property
    def filtered_glass_transcript(self): return FilteredGlassTranscriptBMDC
    @property
    def glass_transcript_source(self): return GlassTranscriptSourceBMDC
    @property
    def glass_transcript_source_prep(self): return GlassTranscriptSourcePrepBMDC
    @property
    def glass_transcript_nucleotides(self): return GlassTranscriptNucleotidesBMDC
    @property
    def glass_transcript_sequence(self): return GlassTranscriptSequenceBMDC
    @property
    def glass_transcript_non_coding(self): return GlassTranscriptNonCodingBMDC
    @property
    def glass_transcript_infrastructure(self): return GlassTranscriptInfrastructureBMDC
    @property
    def glass_transcript_patterned(self): return GlassTranscriptPatternedBMDC
    @property
    def glass_transcript_duped(self): return GlassTranscriptDupedBMDC
    @property
    def glass_transcript_conserved(self): return GlassTranscriptConservedBMDC
    @property
    def glass_transcribed_rna(self): return GlassTranscribedRnaBMDC
    @property
    def glass_transcribed_rna_source(self): return GlassTranscribedRnaSourceBMDC
    @property
    def peak_feature(self): return PeakFeatureBMDC
    @property
    def glass_transcript_label(self): return GlassTranscriptLabelBMDC
    
class GlassTranscriptBMDC(GlassTranscript):
    cell_base = BMDCBase()
    
    #labels = models.ManyToManyField(TranscriptClass, through='GlassTranscriptLabelBMDC')
    
    class Meta:
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Unfiltered Glass transcript (%s)' % CELL_TYPE
        verbose_name_plural = 'Unfiltered Glass transcripts (%s)' % CELL_TYPE

class GlassTranscriptPrepBMDC(GlassTranscriptPrep):
    cell_base = BMDCBase()
    class Meta:
        db_table    = 'glass_atlas_%s_%s_prep"."glass_transcript' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower())
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Preparatory Glass transcript (%s)' % CELL_TYPE
        verbose_name_plural = 'Preparatory Glass transcripts (%s)' % CELL_TYPE

class FilteredGlassTranscriptBMDC(GlassTranscriptBMDC, FilteredGlassTranscript):
    cell_base = BMDCBase()
    objects = FilteredGlassTranscript.objects
    class Meta:
        proxy = True
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name= 'Glass transcript (%s)' % CELL_TYPE
        verbose_name_plural= 'Glass transcripts (%s)' % CELL_TYPE

class GlassTranscriptNucleotidesBMDC(GlassTranscriptNucleotides):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta:
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_nucleotides' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript nucleotide sequence (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript nucleotide sequences (%s)' % CELL_TYPE
      
class GlassTranscriptSourcePrepBMDC(GlassTranscriptSourcePrep):
    glass_transcript = models.ForeignKey(GlassTranscriptPrepBMDC, related_name='glasstranscriptsource')
    cell_base = BMDCBase()
    class Meta:
        db_table    = 'glass_atlas_%s_%s_prep"."glass_transcript_source' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower())
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Preparatory Glass transcript source (%s)' % CELL_TYPE
        verbose_name_plural = 'Preparatory Glass transcript sources (%s)' % CELL_TYPE

class GlassTranscriptSourceBMDC(GlassTranscriptSource):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC, related_name='glasstranscriptsource')
    cell_base = BMDCBase()
    class Meta:
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_source' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript source (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript sources (%s)' % CELL_TYPE
       
class GlassTranscriptSequenceBMDC(GlassTranscriptSequence):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_sequence' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript sequence region (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript sequence regions (%s)' % CELL_TYPE
           
class GlassTranscriptNonCodingBMDC(GlassTranscriptNonCoding):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_non_coding' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript non-coding region (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript non-coding regions (%s)' % CELL_TYPE
        
class GlassTranscriptInfrastructureBMDC(GlassTranscriptInfrastructure):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_infrastructure' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript infrastructure ncRNA region (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript infrastructure ncRNA regions (%s)' % CELL_TYPE
        
class GlassTranscriptPatternedBMDC(GlassTranscriptPatterned):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_patterned' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript patterned region (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript patterned regions (%s)' % CELL_TYPE

class GlassTranscriptDupedBMDC(GlassTranscriptDuped):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_duped' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript duped region (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript duped regions (%s)' % CELL_TYPE
        
class GlassTranscriptConservedBMDC(GlassTranscriptConserved):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_conserved' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript conserved region (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript conserved regions (%s)' % CELL_TYPE

##################################################
# Transcribed RNA
##################################################
class GlassTranscribedRnaBMDC(GlassTranscribedRna):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC, blank=True, null=True)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s_rna%s"."glass_transcribed_rna' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcribed RNA (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcribed RNAs (%s)' % CELL_TYPE
        
class GlassTranscribedRnaSourceBMDC(GlassTranscribedRnaSource):
    glass_transcribed_rna = models.ForeignKey(GlassTranscribedRnaBMDC, related_name='glasstranscribedrnasource')
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s_rna%s"."glass_transcribed_rna_source' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcribed RNA Source (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcribed RNA Sources (%s)' % CELL_TYPE

##################################################
# Features
##################################################       
class PeakFeatureBMDC(PeakFeature):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."peak_feature' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Peak feature (%s)' % CELL_TYPE
        verbose_name_plural = 'Peak feature (%s)' % CELL_TYPE

##################################################
# Features
##################################################       
class GlassTranscriptLabelBMDC(GlassTranscriptLabel):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta:
        db_table    = 'glass_atlas_%s_%s_prep"."glass_transcript_label' % (current_settings.TRANSCRIPT_GENOME, CELL_TYPE.lower())
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript label (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript labels (%s)' % CELL_TYPE