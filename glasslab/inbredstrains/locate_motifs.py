'''
Created on Dec 7, 2011

@author: karmel

Using the genomes generated by generate_genomes_from_referency.py, we can
use Homer to get all known motifs in each of the strains. However, the positions
of the motifs are according to the generated genome, not the reference.

In order to be able to compare motifs across strains, then, we update each motif
that is found by Homer with the appropriate position in the reference genome. 

The logic here goes something like this: The only differences between the reference
and the strain genomes are caused by SNPs/indels. So, for each motif in the strain
genome, we find the SNP/indel that most immediately precedes it. That SNP/indel's
location in the reference genome is already known, so we just count from its
location to where the motif begins to get the motif's start in the reference genome.
'''

from __future__ import division
from glasslab.utils.datatypes.genetics import InbredStrainVariation,\
    InbredStrain
import sys
from glasslab.utils.database import execute_query
from glasslab.utils.datatypes.motif import Motif, BindingSite
from django.db import connection

    
def main(strain='BALB'):
    # First, update the motifs db_table so that we ensure we're looking at the right strain
    motif_db_table = Motif._meta.db_table
    binding_site_db_table = BindingSite._meta.db_table
    Motif._meta.db_table = Motif._meta.db_table.replace('homer_motifs','homer_motifs_%s' % strain.lower()) 
    BindingSite._meta.db_table = BindingSite._meta.db_table.replace('homer_motifs','homer_motifs_%s' % strain.lower())
    
    strain_obj = InbredStrain.objects.get(name=strain)

    for chr_id in xrange(1,23):
        #if chr_id < 21: continue
        print 'Processing chromosome %d.' % chr_id
        
        sites = BindingSite.objects.raw('''
            SELECT "%s_%d"."id", 
            "%s_%d"."chromosome_id", 
            "%s_%d"."motif_id", 
            "%s_%d"."strain_start", 
            "%s_%d"."strain_end", 
            "%s_%d"."strain_start_end"
            FROM "%s_%d" 
            WHERE start IS NULL
            ORDER BY "%s_%d"."strain_start" ASC''' % \
                tuple([BindingSite._meta.db_table, chr_id]*8))
        
        query_string = ''
        last_end, no_variants = None, None # Placeholder for the last variation; convenience so we don't rerun the long query a lot.
        for i,site in enumerate(list(sites)):
            # Get nearest same-strain variant
            var_start = InbredStrainVariation.objects.raw('''
            SELECT "genetics"."inbred_strain_variation_%d"."id", 
            "genetics"."inbred_strain_variation_%d"."chromosome_id", 
            "genetics"."inbred_strain_variation_%d"."strain_start", 
            "genetics"."inbred_strain_variation_%d"."strain_end", 
            "genetics"."inbred_variant_%d"."start" as ref_start
            FROM "genetics"."inbred_strain_variation_%d" 
            INNER JOIN "genetics"."inbred_variant_%d" 
            ON ("genetics"."inbred_strain_variation_%d"."inbred_variant_id" = "genetics"."inbred_variant_%d"."id") 
            WHERE "genetics"."inbred_strain_variation_%d"."inbred_strain_id" = %d
            AND "genetics"."inbred_strain_variation_%d".strain_start <= %d
            ORDER BY "genetics"."inbred_strain_variation_%d"."strain_start" DESC
            LIMIT 1''' % \
                tuple([chr_id]*10 + [strain_obj.id, chr_id, site.strain_start, chr_id]))
            
            # Get nearest end so we can accurately calculate end position as well.
            # Note that we search for the nearest variant end, not start, to cover the case
            # where the variant overlaps with the end of the motif site.
            var_end = InbredStrainVariation.objects.raw('''
            SELECT "genetics"."inbred_strain_variation_%d"."id", 
            "genetics"."inbred_strain_variation_%d"."chromosome_id", 
            "genetics"."inbred_strain_variation_%d"."strain_start", 
            "genetics"."inbred_strain_variation_%d"."strain_end", 
            "genetics"."inbred_variant_%d"."start" as ref_start
            FROM "genetics"."inbred_strain_variation_%d" 
            INNER JOIN "genetics"."inbred_variant_%d" 
            ON ("genetics"."inbred_strain_variation_%d"."inbred_variant_id" = "genetics"."inbred_variant_%d"."id") 
            WHERE "genetics"."inbred_strain_variation_%d"."inbred_strain_id" = %d
            AND "genetics"."inbred_strain_variation_%d".strain_end >= %d
            ORDER BY "genetics"."inbred_strain_variation_%d"."strain_end" ASC
            LIMIT 1''' % \
                tuple([chr_id]*10 + [strain_obj.id, chr_id, site.strain_end, chr_id]))
            
            try:
                var_start = var_start[0]
                var_start_start = var_start.strain_start
                var_start_ref_start = var_start.ref_start
            except IndexError:
                # There are no preceding variations- the alt_start == start in reference
                var_start_start = site.strain_start
                var_start_ref_start = site.strain_start
                
            try:
                var_end = var_end[0]
                var_end_start = var_end.strain_start
                var_end_ref_start = var_end.ref_start
            except IndexError:
                # There are no following variations; get the last available
                
                try:
                    if no_variants: raise IndexError
                    if not last_end:
                        last_end = InbredStrainVariation.objects.raw('''
                        SELECT "genetics"."inbred_strain_variation_%d"."id", 
                        "genetics"."inbred_strain_variation_%d"."chromosome_id", 
                        "genetics"."inbred_strain_variation_%d"."strain_start", 
                        "genetics"."inbred_strain_variation_%d"."strain_end", 
                        "genetics"."inbred_variant_%d"."start" as ref_start
                        FROM "genetics"."inbred_strain_variation_%d" 
                        INNER JOIN "genetics"."inbred_variant_%d" 
                        ON ("genetics"."inbred_strain_variation_%d"."inbred_variant_id" = "genetics"."inbred_variant_%d"."id") 
                        WHERE "genetics"."inbred_strain_variation_%d"."inbred_strain_id" = %d
                        ORDER BY "genetics"."inbred_variant_%d"."end" DESC
                        LIMIT 1''' % \
                            tuple([chr_id]*10 + [strain_obj.id, chr_id]))
                        last_end = last_end[0]
                    var_end = last_end
                    var_end_start = var_end.strain_start
                    var_end_ref_start = var_end.ref_start
                except IndexError:
                    # No variants on this chromosome; alternate is the same as reference
                    no_variants = True
                    var_end_start = site.strain_start
                    var_end_ref_start = site.strain_start
                    
                
            # Update motif record with position in the reference genome
            ref_start = (site.strain_start - var_start_start) + var_start_ref_start
            # Though we used the end to find this variation, we use its start as the adjustment factor
            # to avoid any discrepancies introduced by indels
            ref_end = (site.strain_end - var_end_start) + var_end_ref_start
                
            
            # Run query manually for the sake of speed.
            query_string += '''UPDATE "%s_%d" 
                set "start" = %d,
                    "end" = %d, 
                    "start_end" = public.make_box(%d, 0, %d, 0)
                    WHERE id = %d;
                    ''' % (BindingSite._meta.db_table,
                          chr_id, ref_start, ref_end,
                          ref_start, ref_end, site.id)
            
            # Batch so that we don't have to hit the DB so many times
            if i % 1000 == 0:
                connection.close() 
                execute_query(query_string)
                query_string = ''
                print "Processed motif %d with position %d." % (i, ref_start)
        
        # Get the last batch
        if query_string:
            connection.close() 
            execute_query(query_string)
            query_string = ''
            print "Processed motif %d with position %d." % (i, ref_start)
        
    Motif._meta.db_table = motif_db_table
    BindingSite._meta.db_table = binding_site_db_table
    
if __name__=='__main__':
    if len(sys.argv) > 1: strain = sys.argv[1]
    else: strain = 'BALB'
    print 'Setting reference genome coordinates for motifs for %s.' % strain
    main(strain)
        
                