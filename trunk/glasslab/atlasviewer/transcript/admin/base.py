'''
Created on Nov 8, 2010

@author: karmel
'''
from glasslab.glassatlas.datatypes.transcript import GlassTranscriptSource,\
    GlassTranscriptSequence, GlassTranscriptNonCoding,\
    GlassTranscriptSourcePrep, GlassTranscriptInfrastructure
from glasslab.config import current_settings
from glasslab.atlasviewer.shared.admin import make_all_fields_readonly,\
    ReadOnlyInline, ReadOnlyAdmin
from glasslab.glassatlas.datatypes.feature import PeakFeature

class GlassTranscriptSourceInline(ReadOnlyInline):
    model = GlassTranscriptSource    
    readonly_fields = make_all_fields_readonly(model)

class GlassTranscriptSourcePrepInline(ReadOnlyInline):
    model = GlassTranscriptSourcePrep
    readonly_fields = make_all_fields_readonly(model)

class GlassTranscriptSequenceInline(ReadOnlyInline):
    model = GlassTranscriptSequence
    readonly_fields = make_all_fields_readonly(model)
class GlassTranscriptNonCodingInline(ReadOnlyInline):
    model = GlassTranscriptNonCoding
    readonly_fields = make_all_fields_readonly(model)
class GlassTranscriptInfrastructureInline(ReadOnlyInline):
    model = GlassTranscriptInfrastructure
    readonly_fields = make_all_fields_readonly(model)

class PeakFeatureInline(ReadOnlyInline):
    model = PeakFeature
    readonly_fields = make_all_fields_readonly(model)

class TranscriptBase(ReadOnlyAdmin):
    ordering        = ('chromosome','transcription_start')
    search_fields   = ['transcription_start','transcription_end',]
    
    def transcript_length(self, obj):
        return obj.transcription_end - obj.transcription_start + 1
    transcript_length.short_description = 'Length'
    
    def ucsc_browser_link(self, obj):
        all_tracks_file = 'all_tracks.txt'
        all_link = self._ucsc_browser_link(obj, all_tracks_file, 'All')
        return all_link
                                       
    def _ucsc_browser_link(self, obj, session_file, text): 
                           
        return '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&amp;position=%s%%3A+%d-%d' \
                % (current_settings.GENOME, obj.chromosome.name.strip(), 
                           obj.transcription_start, obj.transcription_end) \
                + '&amp;hgS_doLoadUrl=submit&amp;hgS_loadUrlName=' \
                + obj.ucsc_session_url + session_file \
                + '" target="_blank">View</a>'
                     
    ucsc_browser_link.short_description = 'UCSC Browser' 
    ucsc_browser_link.allow_tags = True 
    
    def truncated_score(self, obj):
        return obj.score is not None and '%.3f' % obj.score or 'None'
    truncated_score.short_description = 'Score' 
    
    def truncated_density(self, obj):
        return obj.density is not None and '%.3f' % obj.density or 'None'
    truncated_density.short_description = 'Density'

class GlassTranscriptAdmin(TranscriptBase):
    list_display    = ('chromosome','transcription_start','transcription_end','strand',
                       'transcript_length', 'truncated_score', 'spliced', 'ucsc_browser_link', 'modified')
    list_filter     = ('chromosome','strand','spliced')
    
    save_on_top     = True
    inlines         = [GlassTranscriptSequenceInline,
                       GlassTranscriptNonCodingInline,
                       GlassTranscriptSourceInline, 
                       ]

class GlassTranscriptPrepAdmin(TranscriptBase):
    list_display    = ('chromosome','transcription_start','transcription_end','strand',
                       'transcript_length', 'ucsc_browser_link')
    list_filter     = ('chromosome','strand')
    inlines         = [GlassTranscriptSourcePrepInline, 
                       ]
        
class PeakFeatureAdmin(ReadOnlyAdmin):
    list_display    = ('glass_transcript_link','relationship','peak_type')
    list_filter     = ('peak_type',)
    
    def glass_transcript_link(self, obj):
        if not obj.glass_transcript: return ''
        return '<a href="/admin/Transcription_%s/glasstranscript%s/%d" target="_blank">%s</a>'\
                            % (obj.cell_base.cell_type, obj.cell_base.cell_type.lower(),
                               obj.glass_transcript.id, str(obj.glass_transcript))
    glass_transcript_link.short_description = 'Glass Transcript' 
    glass_transcript_link.allow_tags = True 
    
    