'''
Created on Sep 27, 2014

@author: karmel
'''
SAMPLES = (('naive', 'atac', ''),
           ('naive', 'atac', 'foxo1_ko_'),
           ('naive', 'foxo1', ''),
           ('lcmv_d12', 'foxo1', ''),
           ('naive', 'h3k4me2', ''),
           ('lcmv_d12', 'h3k4me2', ''))


def sample_name(condition, seq_type, breed):
    return '{}{}_{}'.format(breed, condition, seq_type)


def get_threshold(seq):
    if seq.lower() in ('h3k4me2', ):
        return 20
    else:
        return 10
