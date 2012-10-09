'''
Created on Jan 17, 2011

@author: karmel

Features and annotations relevant to transcripts.
'''
from django.db import models, connection
from glasslab.glassatlas.datatypes.metadata import SequencingRun, PeakType
from glasslab.sequencing.datatypes.tag import wrap_errors
from glasslab.config import current_settings
from glasslab.utils.database import execute_query
from glasslab.glassatlas.datatypes.transcript import multiprocess_all_chromosomes
from glasslab.utils.datatypes.basic_model import GlassModel
from glasslab.sequencing.datatypes.peak import GlassPeak


def wrap_add_features_from_chipseq(cls, chr_list, *args): wrap_errors(cls._add_features_from_chipseq, chr_list, *args)

class PeakFeature(GlassModel):
    # Transcript is cell specific.
    glass_transcript = None
    relationship    = models.CharField(max_length=100)
    peak_type       = models.ForeignKey(PeakType)
    glass_peak_id   = models.IntegerField(max_length=12)
    sequencing_run  = models.ForeignKey(SequencingRun)
    length          = models.IntegerField(max_length=12)
    tag_count       = models.IntegerField(max_length=12)
    distance_to_tss = models.IntegerField(max_length=12)
    
    prep_table = GlassPeak._meta.db_table
    
    class Meta:
        abstract    = True
          
    def __unicode__(self):
        return '{0} {1} {2} (Dist to TSS: {3})'.format(str(self.glass_transcript), 
                             self.relationship.strip(),
                             str(self.peak_type_id),
                             self.distance_to_tss)
    
    
    @classmethod 
    def add_from_peaks(cls,  tag_table):
        connection.close()
        sequencing_run = SequencingRun.objects.get(source_table=tag_table)
        if sequencing_run.type.strip() == 'ChIP-Seq':
            cls.add_features_from_chipseq(tag_table, sequencing_run)
    
    ################################################
    # ChIP-Seq to peak feature
    ################################################
    @classmethod 
    def add_features_from_chipseq(cls,  tag_table, sequencing_run):
        multiprocess_all_chromosomes(wrap_add_features_from_chipseq, cls, sequencing_run)
        
    @classmethod
    def _add_features_from_chipseq(cls, chr_list, sequencing_run):
        for chr_id in chr_list:
            print 'Adding peak features for chromosome %d' % chr_id
            query = """
                SELECT glass_atlas_%s_%s%s.insert_associated_peak_features_from_run(%d, %d);
                """ % (current_settings.GENOME,
                       current_settings.CELL_TYPE.lower(),
                       current_settings.STAGING,
                       sequencing_run.id, chr_id)
            execute_query(query)
    

