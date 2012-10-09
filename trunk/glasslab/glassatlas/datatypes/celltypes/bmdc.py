'''
Created on Dec 22, 2010

@author: karmel
'''
from glasslab.glassatlas.datatypes.transcript import GlassTranscript,\
    FilteredGlassTranscript, GlassTranscriptSource,\
    GlassTranscriptSequence, GlassTranscriptNonCoding, \
    CellTypeBase, GlassTranscriptSourcePrep,\
    GlassTranscriptPrep, GlassTranscriptInfrastructure, \
    FilteredGlassTranscriptManager
from glasslab.config import current_settings
from django.db import models
from glasslab.glassatlas.datatypes.feature import PeakFeature

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
    def glass_transcript_sequence(self): return GlassTranscriptSequenceBMDC
    @property
    def glass_transcript_non_coding(self): return GlassTranscriptNonCodingBMDC
    @property
    def glass_transcript_infrastructure(self): return GlassTranscriptInfrastructureBMDC
    @property
    def peak_feature(self): return PeakFeatureBMDC
    
class GlassTranscriptBMDC(GlassTranscript):
    cell_base = BMDCBase()
    
    #labels = models.ManyToManyField(TranscriptClass, through='GlassTranscriptLabelBMDC')
    
    class Meta:
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript' % (current_settings.GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Unfiltered Glass transcript (%s)' % CELL_TYPE
        verbose_name_plural = 'Unfiltered Glass transcripts (%s)' % CELL_TYPE

class GlassTranscriptPrepBMDC(GlassTranscriptPrep):
    cell_base = BMDCBase()
    class Meta:
        db_table    = 'glass_atlas_%s_%s_prep"."glass_transcript' % (current_settings.GENOME, CELL_TYPE.lower())
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Preparatory Glass transcript (%s)' % CELL_TYPE
        verbose_name_plural = 'Preparatory Glass transcripts (%s)' % CELL_TYPE

class FilteredGlassTranscriptBMDC(GlassTranscriptBMDC, FilteredGlassTranscript):
    cell_base = BMDCBase()
    objects = FilteredGlassTranscriptManager()
    class Meta:
        proxy = True
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name= 'Glass transcript (%s)' % CELL_TYPE
        verbose_name_plural= 'Glass transcripts (%s)' % CELL_TYPE

class GlassTranscriptSourcePrepBMDC(GlassTranscriptSourcePrep):
    glass_transcript = models.ForeignKey(GlassTranscriptPrepBMDC, related_name='glasstranscriptsource')
    cell_base = BMDCBase()
    class Meta:
        db_table    = 'glass_atlas_%s_%s_prep"."glass_transcript_source' % (current_settings.GENOME, CELL_TYPE.lower())
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Preparatory Glass transcript source (%s)' % CELL_TYPE
        verbose_name_plural = 'Preparatory Glass transcript sources (%s)' % CELL_TYPE

class GlassTranscriptSourceBMDC(GlassTranscriptSource):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC, related_name='glasstranscriptsource')
    cell_base = BMDCBase()
    class Meta:
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_source' % (current_settings.GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript source (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript sources (%s)' % CELL_TYPE
       
class GlassTranscriptSequenceBMDC(GlassTranscriptSequence):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_sequence' % (current_settings.GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript sequence region (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript sequence regions (%s)' % CELL_TYPE
           
class GlassTranscriptNonCodingBMDC(GlassTranscriptNonCoding):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_non_coding' % (current_settings.GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript non-coding region (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript non-coding regions (%s)' % CELL_TYPE
        
class GlassTranscriptInfrastructureBMDC(GlassTranscriptInfrastructure):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."glass_transcript_infrastructure' % (current_settings.GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Glass transcript infrastructure ncRNA region (%s)' % CELL_TYPE
        verbose_name_plural = 'Glass transcript infrastructure ncRNA regions (%s)' % CELL_TYPE
        
##################################################
# Features
##################################################       
class PeakFeatureBMDC(PeakFeature):
    glass_transcript = models.ForeignKey(GlassTranscriptBMDC)
    cell_base = BMDCBase()
    class Meta: 
        db_table    = 'glass_atlas_%s_%s%s"."peak_feature' % (current_settings.GENOME, CELL_TYPE.lower(), current_settings.STAGING)
        app_label   = 'Transcription_%s' % CELL_TYPE
        verbose_name = 'Peak feature (%s)' % CELL_TYPE
        verbose_name_plural = 'Peak feature (%s)' % CELL_TYPE

