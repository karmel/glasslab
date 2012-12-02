'''
Created on Sep 24, 2010

@author: karmel

This module serves as a singleton settings object, for settings that
should be set on a run-wide basis.
'''
#####################################
# Genomes
#####################################
GENOME_CHOICES = ['mm9']

GENOME = 'mm9'
GENOME_CHROMOSOMES = range(1,23)


CELL_TYPE = 'ThioMac'

STAGING = '' # Set to the appropriate suffix during DB staging.
STAGING_SUFFIX = '_staging'

GENOME_ASSEMBLY_PATHS = {'mm9': '/Volumes/Unknowme/kallison/Genomes/mm9/fasta',}

#####################################
# Databases
#####################################
CURRENT_SCHEMA = 'current_projects'

#####################################
# Compute power
#####################################
ALLOWED_PROCESSES = 6
CHR_LISTS = None # Dynamically set during processing

#####################################
# Compute resources
#####################################
# Used to log in as postgres user. Note that authorized keys must be set up!
PG_ACCESS_CMD = 'ssh postgres@glass.bioinforma.tc' 
PG_HOME = '/Library/PostgreSQL/9.1'
