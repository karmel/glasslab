'''
Created on Dec 22, 2010

@author: karmel
'''
from glasslab.atlasviewer.transcript.admin.base import GlassTranscriptAdmin,\
    GlassTranscriptSequenceInline,\
    GlassTranscriptNonCodingInline, GlassTranscriptSourceInline,\
    PeakFeatureInline, PeakFeatureAdmin,\
    GlassTranscriptPrepAdmin, GlassTranscriptSourcePrepInline,\
    GlassTranscriptInfrastructureInline
from glasslab.glassatlas.datatypes.celltypes.thiomac import FilteredGlassTranscriptThioMac,\
    GlassTranscriptThioMac, GlassTranscriptInfrastructureThioMac,\
    GlassTranscriptSequenceThioMac, GlassTranscriptSourceThioMac,\
    GlassTranscriptNonCodingThioMac, PeakFeatureThioMac,\
    GlassTranscriptSourcePrepThioMac, GlassTranscriptPrepThioMac
from django.contrib import admin


class GlassTranscriptSequenceThioMacInline(GlassTranscriptSequenceInline):
    model = GlassTranscriptSequenceThioMac
class GlassTranscriptNonCodingThioMacInline(GlassTranscriptNonCodingInline):
    model = GlassTranscriptNonCodingThioMac
class GlassTranscriptInfrastructureThioMacInline(GlassTranscriptInfrastructureInline):
    model = GlassTranscriptInfrastructureThioMac
class GlassTranscriptSourceThioMacInline(GlassTranscriptSourceInline):
    model = GlassTranscriptSourceThioMac

class GlassTranscriptSourcePrepThioMacInline(GlassTranscriptSourcePrepInline):
    model = GlassTranscriptSourcePrepThioMac

class PeakFeatureThioMacInline(PeakFeatureInline):
    model = PeakFeatureThioMac

class GlassTranscriptThioMacAdmin(GlassTranscriptAdmin):
    search_fields   = ['transcription_start','transcription_end']
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

