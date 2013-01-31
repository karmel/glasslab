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
from glasslab.dataanalysis.graphing.seq_grapher import SeqGrapher
import numpy

kla_col='kla_6h_lfc'
    
if __name__ == '__main__':
    yzer = SeqGrapher()
    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/HiC/'
    dirpath = yzer.get_path(dirpath)
    data_dirpath = yzer.get_filename(dirpath, 'enhancer_sets')
    img_dirpath = yzer.get_and_create_path(dirpath, 'genes_to_average_enhancer_lfc')
    
    keys = ('all', 'notx', 'kla', 'notx_only', 'kla_only', 'shared_enh')
    
    if True:
        interactions = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_pairs_refseq_with_me2.txt'))
        interactions = interactions[interactions['count'] > 1]
    
        all_transcripts = yzer.import_file(yzer.get_filename(data_dirpath,'transcript_vectors.txt'))
        
        
        transcripts = all_transcripts[['id', 'kla_lfc',  'kla_6h_lfc']]
            
        # Associate gene id
        interactions = interactions.merge(transcripts, how='left', on='id')
        
        transcripts['id_2'] = transcripts['id']
        transcripts = transcripts.drop(['id'], axis=1)
        interactions = interactions.merge(transcripts, how='left', on='id_2', suffixes=['','_2'])
        
        interactions = interactions.fillna(0)
        
        # Key on peak id, not enhancer id, which could be bidirectional
        #interactions['id_2'] = interactions['h3k4me2_id']
        interactions['hash'] = interactions.apply(lambda row: '{0}.{1}'.format(row['id'],row['id_2']), axis=1)
        
        
        # Group by gene
        data = interactions.groupby(by='id')
        
        keys = ('all', 'notx', 'kla', 'notx_only', 'kla_only', 'shared_enh')
        expression_vectors = []
        for gene_id, group in data:
            notx = group[group['sequencing_run_id'] == 765].sort(['norm_count', kla_col], ascending=False)
            kla_30m = group[group['sequencing_run_id'] == 766].sort(['norm_count', kla_col], ascending=False)
            kla_4h = group[group['sequencing_run_id'] == 773].sort(['norm_count', kla_col], ascending=False)
            
            if len(group) == len(kla_30m): continue # KLA 30m only
            
            shared = kla_4h[(kla_4h['hash'].isin(notx['hash']))]['hash']
            new_ixns = group[~group['hash'].isin(shared)]
            
            all_enh = group[kla_col + '_2'].mean()
            notx = notx[kla_col + '_2'].mean()
            kla = kla_4h[kla_col + '_2'].mean()
            notx_only = new_ixns[new_ixns['sequencing_run_id'] == 765][kla_col + '_2'].mean()
            kla_only = new_ixns[new_ixns['sequencing_run_id'] == 773][kla_col + '_2'].mean()
            shared_enh = group[group['hash'].isin(shared)][kla_col + '_2'].mean()
            
            expression_vectors.append([gene_id, group[kla_col].mean(), 
                                       all_enh, notx, kla, notx_only, kla_only, shared_enh])
            
        
        # Sort by gene LFD
        expression_vectors.sort(key=lambda row: row[1], reverse=True)
        
        # Print out as cdt file, as if it were clustered, for feeding into JavaTree
        # Note that we print id, weight, gene, enh
        for i, name in enumerate(keys):
            outfile = open(yzer.get_filename(img_dirpath, name + '_vectors.cdt'), 'w')
            outfile.write('id\tweight\tgene_lfc\tmean_enhancer_lfc\n')
            for line in expression_vectors: 
                if not numpy.isnan(line[i+2]):
                    outfile.write('{0}\t1.0\t{1}\t{2}\n'.format(line[0], line[1], line[i+2]))
            outfile.close()
    
    for i, name in enumerate(keys):
        data = yzer.import_file(yzer.get_filename(img_dirpath, name + '_vectors.cdt'))
        data = data[data['gene_lfc'].abs() > 1]
        print name, len(data)
        print data.corr()
        