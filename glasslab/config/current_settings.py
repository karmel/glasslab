'''
Created on Sep 24, 2010

@author: karmel

This module serves as a singleton settings object, for settings that
should be set on a run-wide basis.
'''

#####################################
# Genomes
#####################################
GENOME_CHOICES = {'mm9': {'name':'Mus musculus', 'chromosomes': range(1,23)},
                  'dm3': {'name':'Drosophila melanogaster', 'chromosomes': range(1,15)}}

GENOME = 'mm9'


CELL_TYPE = 'Default'

STAGING = '' # Set to the appropriate suffix during DB staging.
STAGING_SUFFIX = '_staging'

MAX_EDGE = 100 # Max edge length in 2D between two proto-transcripts

UCSC_SESSION = 'TCC_2012'
#####################################
# Databases
#####################################
CURRENT_SCHEMA = 'current_projects'
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'glasslab',
        'USER': 'glass',
        'PASSWORD': 'monocyte',
        'HOST': 'localhost',
        'PORT': '54321',
    },
    'read_only': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'glasslab',
        'USER': 'glass_read_only',
        'PASSWORD': 'monocyte',
        'HOST': 'localhost',
        'PORT': '54321',
    }
}

#####################################
# Compute power
#####################################
ALLOWED_PROCESSES = 6
CHR_LISTS = None  # Dynamically set during processing

#####################################
# Compute resources
#####################################
# Used to log in as postgres user. Note that authorized keys must be set up!
PG_ACCESS_CMD = 'ssh postgres@glass.bioinforma.tc' 
PG_HOME = '/Library/PostgreSQL/9.1'

#####################################
# Required for Django; not used
#####################################
SECRET_KEY = 'feg^reh@(rdyue(yfawu0mg532ok^yfl9$1%*ge+ng$1@0gf%x'
