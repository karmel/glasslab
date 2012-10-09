'''
Created on Dec 22, 2010

@author: karmel
'''
from glasslab.atlasviewer.transcript.admin.base import GlassTranscriptAdmin,\
    GlassTranscriptSequenceInline,\
    GlassTranscriptNonCodingInline, GlassTranscriptSourceInline,\
    PeakFeatureInline, PeakFeatureAdmin,\
    GlassTranscriptPrepAdmin, GlassTranscriptSourcePrepInline,\
    GlassTranscriptDupedInline, GlassTranscriptInfrastructureInline
from glasslab.glassatlas.datatypes.celltypes.thiomac import FilteredGlassTranscriptThioMac,\
    GlassTranscriptThioMac, GlassTranscriptDupedThioMac, \
    GlassTranscriptSequenceThioMac, GlassTranscriptSourceThioMac,\
    GlassTranscriptNonCodingThioMac, PeakFeatureThioMac,\
    GlassTranscriptSourcePrepThioMac, GlassTranscriptPrepThioMac,\
    GlassTranscriptInfrastructureThioMac
from django.contrib import admin


class GlassTranscriptSequenceThioMacInline(GlassTranscriptSequenceInline):
    model = GlassTranscriptSequenceThioMac
class GlassTranscriptNonCodingThioMacInline(GlassTranscriptNonCodingInline):
    model = GlassTranscriptNonCodingThioMac
class GlassTranscriptDupedThioMacInline(GlassTranscriptDupedInline):
    model = GlassTranscriptDupedThioMac
class GlassTranscriptInfrastructureThioMacInline(GlassTranscriptInfrastructureInline):
    model = GlassTranscriptInfrastructureThioMac
class GlassTranscriptSourceThioMacInline(GlassTranscriptSourceInline):
    model = GlassTranscriptSourceThioMac

class GlassTranscriptSourcePrepThioMacInline(GlassTranscriptSourcePrepInline):
    model = GlassTranscriptSourcePrepThioMac

class PeakFeatureThioMacInline(PeakFeatureInline):
    model = PeakFeatureThioMac

class GlassTranscriptThioMacAdmin(GlassTranscriptAdmin):
    search_fields   = ['glasstranscriptlabelthiomac__transcript_class__label','transcription_start','transcription_end']
    inlines         = [GlassTranscriptSequenceThioMacInline,
                       GlassTranscriptNonCodingThioMacInline,
                       GlassTranscriptInfrastructureThioMacInline,
                       GlassTranscriptSourceThioMacInline, 
                       PeakFeatureThioMacInline, 
                       ]
    
class GlassTranscriptPrepThioMacAdmin(GlassTranscriptPrepAdmin):
    inlines         = [GlassTranscriptSourcePrepThioMacInline, 
                       ]

class PeakFeatureThioMacAdmin(PeakFeatureAdmin):
    pass

admin.site.register(FilteredGlassTranscriptThioMac, GlassTranscriptThioMacAdmin)
admin.site.register(GlassTranscriptThioMac, GlassTranscriptThioMacAdmin)
admin.site.register(GlassTranscriptPrepThioMac, GlassTranscriptPrepThioMacAdmin)
admin.site.register(PeakFeatureThioMac, PeakFeatureThioMacAdmin)

