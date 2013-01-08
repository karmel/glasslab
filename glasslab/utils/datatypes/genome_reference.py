'''
Created on Sep 24, 2010

@author: karmel

.. attention::
    
    All sequence indices assume a start position of 0, as per the UCSC Genome Browser standard.

'''
from django.db import models
from glasslab.utils.datatypes.basic_model import GlassModel
    
#######################################################
# Genome identifiers
#######################################################

class GenomeType(GlassModel):
    '''
    Genome types with unique gene records.
    '''
    genome_type = models.CharField(max_length=20)
    
    class Meta: 
        db_table    = 'genome_reference"."genome_type'
        app_label   = 'Genome_Reference'

class Genome(GlassModel):
    '''
    Genomes for which we store data.
    '''
    genome_type = models.ForeignKey(GenomeType)
    genome      = models.CharField(max_length=10)
    description = models.CharField(max_length=50)
    
    class Meta: 
        db_table    = 'genome_reference"."genome'
        app_label   = 'Genome_Reference'

class KeggPathway(GlassModel):
    '''
    Kegg Pathway descriptions
    '''
    pathway_key = models.CharField(max_length=10, help_text='Identifier string for the pathway (i.e., "mmu00051"')
    description = models.CharField(max_length=255)
    
    class Meta: 
        db_table    = 'genome_reference"."kegg_pathway'
        app_label   = 'Genome_Reference'
    
    def __unicode__(self): return self.description
    
    

