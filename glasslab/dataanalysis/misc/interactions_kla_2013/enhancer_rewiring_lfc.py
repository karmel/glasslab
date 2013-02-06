'''
Created on Jan 30, 2013

@author: karmel

We want to generate sets of vectors for LFC correlation.
Each vector should have the LFC of:
    enhancer    gene assoc in notx    gene assoc in KLA

So we group by enhancer, and then within each enhancer group,
take the kla_6h_lfc val for the notx gene and the KLA gene.

The complications are enhancers with more than one gene in each
condition, or with an unequal number of genes in each condition.
For those, we will sort genes in each condition by number
of interactions, and allow for null values when there is a number
mismatch.
'''
from __future__ import division
from pandas import DataFrame, Series
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher

kla_col='p65_tag_count'
    
def get_vectors_for_group(group, f_condition=None):
    notx = group[group['sequencing_run_id'] == 765].sort(['norm_count', kla_col], ascending=False)
    kla_4h = group[group['sequencing_run_id'] == 773].sort(['norm_count', kla_col], ascending=False)
    
    df = DataFrame({'enhancer_id': Series(),
                    'enhancer_lfc': Series(),
                    'kla_4h_lfc': Series(kla_4h[kla_col], index=range(len(kla_4h))),
                    'notx_lfc': Series(notx[kla_col], index=range(len(notx))),
                    },)
    df['enhancer_id'] = group['id_2'].mean()
    df['enhancer_lfc'] = group['p65_tag_count_2'].mean()
    if f_condition: df = df[f_condition(df)]
    return df

if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_sets')
    img_dirpath = yzer.get_and_create_path(dirpath, 'enhancer_rewiring_lfc','p65_tags')
    
    interactions = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_refseq_with_me2.txt'))
    interactions = interactions[interactions['count'] > 1]

    transcripts = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_vectors.txt'))
    
    transcripts['kla_6h_rpbp'] = transcripts['kla_6h_tag_count']/(transcripts['length'])*1000
    transcripts['kla_rpbp'] = transcripts['kla_tag_count']/(transcripts['length'])*1000
        
    # Associate gene id
    interactions = interactions.merge(transcripts, how='left', on='id')
    
    transcripts['id_2'] = transcripts['id']
    transcripts = transcripts.drop(['id'], axis=1)
    interactions = interactions.merge(transcripts, how='left', on='id_2', suffixes=['','_2'])
    
    interactions = interactions.fillna(0)
    
    # Key on peak id, not enhancer id, which could be bidirectional
    #interactions['id_2'] = interactions['h3k4me2_id']
    interactions['hash'] = interactions.apply(lambda row: '{0}.{1}'.format(row['id'],row['id_2']), axis=1)
    
    
    # Group by enhancer
    data = interactions.groupby(by='id_2')
    
    all_vectors, rewired, notx_only, kla_only, shared_vectors = [], [], [], [], []
    for enh_id, group in data:
        notx = group[group['sequencing_run_id'] == 765].sort(['norm_count', kla_col], ascending=False)
        kla_30m = group[group['sequencing_run_id'] == 766].sort(['norm_count', kla_col], ascending=False)
        kla_4h = group[group['sequencing_run_id'] == 773].sort(['norm_count', kla_col], ascending=False)
        
        if len(group) == len(kla_30m): continue # KLA 30m only
        
        df = get_vectors_for_group(group)
        
        all_vectors += df.as_matrix().tolist()
    
        # Now do only enhancer pairings that are not shared.
        shared = kla_4h[(kla_4h['hash'].isin(notx['hash']))]['hash']
        new_ixns = group[~group['hash'].isin(shared)]
        
        # All rewired
        df = get_vectors_for_group(new_ixns)
        if len(df) > 0:
            rewired += df.as_matrix().tolist()
        
        # Only those with notx intxns
        df = get_vectors_for_group(new_ixns, lambda g: ~g['notx_lfc'].isnull())
        if len(df) > 0:
            notx_only += df.as_matrix().tolist()
        # Only those with kla intxns
        df = get_vectors_for_group(new_ixns, lambda g: ~g['kla_4h_lfc'].isnull())
        if len(df) > 0:
            kla_only += df.as_matrix().tolist()
        
        # Now do only enhancer pairings that are shared.
        df = get_vectors_for_group(group[group['hash'].isin(shared)])
        if len(df) > 0: 
            shared_vectors += df.as_matrix().tolist()
        
    for name, vectors in (('all',all_vectors),
                          ('rewired', rewired),
                          ('notx_only', notx_only),
                          ('kla_only', kla_only),
                          ('shared', shared_vectors)
                          ):
        # Sort by enhancer LFD
        vectors.sort(key=lambda row: row[1], reverse=True)
        
        # Print out as cdt file, as if it were clustered, for feeding into JavaTree
        # Note that we print id, weight, enh, notx, kla; that requires swapping the 
        # columns that pandas alphabetizes
        outfile = open(yzer.get_filename(img_dirpath, name + '_vectors.cdt'), 'w')
        outfile.write('id\tweight\tenhancer_lfc\tnotx_gene_lfc\tkla_4h_gene_lfc\n')
        for line in vectors: outfile.write('{0}\t1.0\t{1}\t{3}\t{2}\n'.format(*line))
    