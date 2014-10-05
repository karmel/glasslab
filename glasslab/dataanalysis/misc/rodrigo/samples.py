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


# Dataset 2
def get_breed_sets():
    # Set up names
    timepoints = (0, 7, 12, 35)
    celltypes = ('klrghi', 'klrglo')
    breeds = ('', 'foxo1_ko_')

    breed_sets = []
    for breed in breeds:
        samples = []
        short_names = []

        for tp in timepoints:
            if tp == 0:
                samples.append('{}cd8tcell_atac_lcmv_{}d'.format(
                    breed, tp))
                short_names.append('{}d{}'.format(
                    breed, tp))
            else:
                for ct in celltypes:
                    samples.append('{}cd8tcell_{}_atac_lcmv_{}d'.format(
                        breed, ct, tp))
                    short_names.append('{}{}_d{}'.format(
                        breed, ct, tp))
        breed_sets.append([samples, short_names])

    return breed_sets
