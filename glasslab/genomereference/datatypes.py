'''
Created on Sep 24, 2010

@author: karmel

.. attention::
    
    All sequence indices assume a start position of 0, as per the UCSC Genome Browser standard.

'''
from django.db import models
from glasslab.utils.datatypes.basic_model import BoxField, GlassModel
from glasslab.utils.datatypes.genome_reference import KeggPathway
from glasslab.config import current_settings

SCHEMA_BASE = 'genome_reference_{0}'.format(current_settings.GENOME)


class GenomeReferenceBase(GlassModel):
    schema_base = SCHEMA_BASE
    class Meta:
        abstract = True
        
#######################################################
# Per-genome Gene identifiers
#######################################################
class Chromosome(GenomeReferenceBase):
    '''
    Unique record of chromosome, i.e. 'chr1', 'chrUn_random', etc
    '''
    name = models.CharField(max_length=25, blank=False)
    length = models.IntegerField(max_length=25, blank=False)
    
    class Meta:
        db_table = '{0}"."chromosome'.format(SCHEMA_BASE)
        app_label   = 'Genome_Reference'
    
        
    def __unicode__(self): return self.name

class SequenceIdentifier(GenomeReferenceBase):
    '''
    Gene and sequence (i.e., noncoding RNA) identifiers from RefSeq, unique per genome.
    '''
    sequence_identifier = models.CharField(max_length=50, blank=False)
    
    class Meta: 
        db_table = '{0}"."sequence_identifier'.format(SCHEMA_BASE)
        app_label   = 'Genome_Reference'
    
    def __unicode__(self): return self.sequence_identifier
    
    _sequence_detail = None
    @property 
    def sequence_detail(self):
        if not self._sequence_detail:
            detail = SequenceDetail.objects.filter(sequence_identifier=self).order_by('-gene_name')[:1]
            if detail: self._sequence_detail =  detail[0]
        return self._sequence_detail
    
    _sequence_transcription_region = None
    @property 
    def sequence_transcription_region(self):
        if not self._sequence_transcription_region:
            reg = SequenceTranscriptionRegion.objects.filter(sequence_identifier=self).order_by('transcription_start')[:1]
            if reg: self._sequence_transcription_region =  reg[0]
        return self._sequence_transcription_region

class SequenceDetail(GenomeReferenceBase):
    '''
    Gene details, keyed to unique sequences.
    '''
    sequence_identifier = models.ForeignKey(SequenceIdentifier)
    gene_name           = models.CharField(max_length=100, blank=True)
    description         = models.CharField(max_length=255, blank=True)
    ensembl_id          = models.CharField(max_length=100, blank=True)
    pfam_id             = models.CharField(max_length=100, blank=True)
    
    class Meta: 
        db_table    = '{0}"."sequence_detail'.format(SCHEMA_BASE)
        app_label   = 'Genome_Reference'
    
    def __unicode__(self): 
        return '{0} ({1})'.format(self.gene_name, self.sequence_identifier.sequence_identifier)

class SequenceKeggPathway(GlassModel):
    '''
    Mappings of transcription regions and coding sites.
    '''
    sequence_identifier = models.ForeignKey(SequenceIdentifier)
    kegg_pathway        = models.ForeignKey(KeggPathway)
    map_location        = models.CharField(max_length=50, help_text='Mappable identifier for this sequence and pathway; can be used in Kegg URLs.')
    
    class Meta: 
        db_table    = '{0}"."sequence_kegg_pathway'.format(SCHEMA_BASE)
        app_label   = 'Genome_Reference'

class NonCodingRna(GenomeReferenceBase):
    '''
    Unique name and type for ncRNA
    
    '''
    type                = models.CharField(max_length=20)
    description         = models.CharField(max_length=100)
    
    class Meta: 
        db_table    = '{0}"."non_coding_rna'.format(SCHEMA_BASE)
        app_label   = 'Genome_Reference'
        verbose_name = 'Non coding RNA'
        
    def __unicode__(self):
        return '{0} {1}'.format(self.type, self.description.strip())
    
    _non_coding_transcription_region = None
    @property 
    def non_coding_transcription_region(self):
        if not self._non_coding_transcription_region:
            reg = NonCodingTranscriptionRegion.objects.filter(non_coding_rna=self).order_by('transcription_start')[:1]
            if reg: self._non_coding_transcription_region =  reg[0]
        return self._non_coding_transcription_region

#######################################################
# Chromosome region details 
#######################################################
class SequenceTranscriptionRegion(GenomeReferenceBase):
    '''
    Mappings of transcription regions and coding sites.
    '''
    sequence_identifier = models.ForeignKey(SequenceIdentifier)
    chromosome          = models.ForeignKey(Chromosome)
    bin                 = models.IntegerField(max_length=5, help_text='Base-2 determined bin.')
    strand              = models.IntegerField(max_length=1, help_text='0 for +, 1 for -')
    transcription_start = models.IntegerField(max_length=12)
    transcription_end   = models.IntegerField(max_length=12)    
    coding_start        = models.IntegerField(max_length=12)
    coding_end          = models.IntegerField(max_length=12)
    
    start_end           = BoxField(null=True, default=None, help_text='This is a placeholder for the PostgreSQL box type.')
    
    table_name = 'sequence_transcription_region'
    class Meta: 
        db_table    = '{0}"."sequence_transcription_region'.format(SCHEMA_BASE)
        app_label   = 'Genome_Reference'

    def __unicode__(self):
        return 'Sequence Transcription Region for {0}'.format(self.sequence_identifier.sequence_identifier.strip())
        
class SequenceExon(GlassModel):
    '''
    Mappings of transcription regions and coding sites.
    '''
    sequence_transcription_region = models.ForeignKey(SequenceTranscriptionRegion)
    exon_start = models.IntegerField(max_length=12)
    exon_end   = models.IntegerField(max_length=12)    
    frame      = models.IntegerField(max_length=5, help_text='Number o nucleotides needed from prior exon to make a complete amino acid at the start of this exon.')
    start_end  = BoxField(null=True, default=None, help_text='This is a placeholder for the PostgreSQL box type.')
    class Meta: 
        db_table    = '{0}"."sequence_exon'.format(SCHEMA_BASE)
        app_label   = 'Genome_Reference'
        
class NonCodingTranscriptionRegion(GenomeReferenceBase):
    '''
    Mappings of transcription regions that are not tied to RefSeq genes.
    
    Scores are 0 - 1000, higher indicating that the sequence is more likely to be true ncRNA.
    
    '''
    non_coding_rna      = models.ForeignKey(NonCodingRna)
    chromosome          = models.ForeignKey(Chromosome)
    strand              = models.IntegerField(max_length=1, help_text='0 for +, 1 for -')
    transcription_start = models.IntegerField(max_length=12)
    transcription_end   = models.IntegerField(max_length=12)    
    score               = models.IntegerField(max_length=5)
    
    start_end           = BoxField(null=True, default=None, help_text='This is a placeholder for the PostgreSQL box type.')
    
    class Meta: 
        db_table    = '{0}"."non_coding_transcription_region'.format(SCHEMA_BASE)
        app_label   = 'Genome_Reference'

    def __unicode__(self):
        return '{0} Transcription Region for {1}'.format(self.non_coding_rna.type.strip(), self.non_coding_rna.description.strip())

class InfrastructureTranscriptionRegion(GlassModel):
    '''
    Mappings of patterns-- i.e., repeats-- onto transcription regions.
    '''
    type                = models.CharField(max_length=20)
    name                = models.CharField(max_length=100)
    chromosome          = models.ForeignKey(Chromosome)
    strand              = models.IntegerField(max_length=1, help_text='0 for +, 1 for -. Default NULL')
    transcription_start = models.IntegerField(max_length=12)
    transcription_end   = models.IntegerField(max_length=12)
    
    start_end           = BoxField(null=True, default=None, help_text='This is a placeholder for the PostgreSQL box type.')
    
    class Meta: 
        db_table    = '{0}"."patterned_transcription_region'.format(SCHEMA_BASE)
        app_label   = 'Genome_Reference'

    def __unicode__(self):
        return 'Infrastructure Transcription Region for %s %s' % (self.type, self.name.strip())


#######################################################
# Metadata
####################################################### 
ALT_SCHEMA_BASE = 'glass_atlas_mm9'  
class PeakType(GlassModel):
    type    = models.CharField(max_length=50)
    diffuse = models.BooleanField(default=False)
    
    class Meta:
        db_table    = '{0}"."peak_type'.format(ALT_SCHEMA_BASE)
        app_label   = 'Genome_Reference' 
          
    def __unicode__(self):
        return self.type.strip()
    
class SequencingRun(GlassModel):
    '''
    Record of details of a given sequencing run and its total tags.
    '''
    type            = models.CharField(max_length=50, choices=[(x,x) for x in ('Gro-Seq','RNA-Seq','ChIP-Seq','Ribo-Seq','HiC')], 
                                       default='Gro-Seq')
    
    cell_type       = models.CharField(max_length=50)
    name            = models.CharField(max_length=100)
    source_table    = models.CharField(max_length=100)
    description     = models.CharField(max_length=255, blank=True)
    total_tags      = models.IntegerField(max_length=12, help_text='Or peaks if this is a PeakType table.')
    percent_mapped  = models.DecimalField(max_digits=5, decimal_places=2, blank=True, null=True, default=None,
                                          help_text='What percent of tags were successfully mapped by Bowtie?')
    
    peak_type       = models.ForeignKey(PeakType, null=True, default=None, blank=True, help_text='Does this run produce ChIP peaks?')
    
    standard        = models.BooleanField(default=False, 
                                          help_text='Is this a standard run, that should affect expected scores and transcript boundaries?')
    use_for_scoring = models.BooleanField(default=False, 
                                          help_text='Should this run be considered when scoring transcripts?')
    requires_reload = models.BooleanField(default=False, help_text='Should this run be reloaded with the next processing cycle?')
    
    timepoint       = models.CharField(max_length=10, blank=True, null=True)
    wt              = models.BooleanField(default=False)    
    notx            = models.BooleanField(default=False)    
    kla             = models.BooleanField(default=False)    
    other_conditions= models.BooleanField(default=False)  
    
    strain          = models.CharField(max_length=100, blank=True, null=True, help_text='Mouse strain')
    modified        = models.DateTimeField(auto_now=True)
    created         = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        db_table = '{0}"."sequencing_run'.format(ALT_SCHEMA_BASE)
        app_label   = 'Genome_Reference'
        
    def __unicode__(self):
        return '{0} ({1}, "{2}")'.format(self.name, self.type, self.source_table.strip())
        
class SequencingRunAnnotation(GlassModel):
    '''
    Various freeform notes that can be attached to sequencing runs.
    '''
    sequencing_run  = models.ForeignKey(SequencingRun)
    note            = models.CharField(max_length=100)
     
    class Meta:
        db_table    = '{0}"."sequencing_run_annotation'.format(ALT_SCHEMA_BASE)
        app_label   = 'Genome_Reference'
        
    def __unicode__(self):
        return '"{0}" note: {1}'.format(self.sequencing_run.source_table.strip(), self.note)