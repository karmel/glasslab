'''
Created on Nov 16, 2011

@author: karmel
'''

from glasslab.utils.datatypes.basic_model import GlassModel, BoxField
from django.db import models
from glasslab.utils.datatypes.genome_reference import Chromosome

class InbredStrain(GlassModel):
    name = models.CharField(max_length=10)
    long_name = models.CharField(max_length=10)
    
    class Meta: 
        db_table    = 'genetics"."inbred_strain'
        app_label   = 'Genetics'


class InbredVariant(GlassModel):
    type = models.CharField(max_length=10, choices=[(x,x) for x in ('SNP','insertion','deletion','CNV')])
    chromosome = models.ForeignKey(Chromosome)
    start = models.IntegerField(max_length=12)
    end = models.IntegerField(max_length=12)
    reference = models.CharField(max_length=255)
    
    start_end = BoxField(null=True, default=None, help_text='This is a placeholder for the PostgreSQL box type.')
    
    class Meta: 
        db_table    = 'genetics"."inbred_variant'
        app_label   = 'Genetics'

class InbredStrainVariation(GlassModel):
    chromosome = models.ForeignKey(Chromosome)
    inbred_strain = models.ForeignKey(Chromosome)
    inbred_variant = models.ForeignKey(Chromosome)
    alternate = models.CharField(max_length=255)
    
    class Meta: 
        db_table    = 'genetics"."inbred_strain_variation'
        app_label   = 'Genetics'
