'''
Created on Nov 7, 2012

@author: karmel
'''
comparison_sets = (('wt_nc1', 'wt_hfd1', 1.087672),
                    ('wt_nc3', 'wt_hfd3', 0.935613),
                    
                    ('ldlr_ko_nc4', 'ldlr_ko_hfd4', 48.267971),
                    ('ldlr_ko_nc6', 'ldlr_ko_hfd6', 3.116647),
                    
                    ('igfii_ldlr_ko_nc5', 'igfii_ldlr_ko_hfd5', 0.984142),
                    ('igfii_ldlr_ko_nc7', 'igfii_ldlr_ko_hfd7', 1.199993),
                    
                    ('wt_nc1', 'wt_nc3', 2.042524),
                    ('wt_hfd1', 'wt_hfd3', 1.665349),
                    
                    ('ldlr_ko_nc4', 'ldlr_ko_nc6', 0.952442),
                    ('ldlr_ko_hfd4', 'ldlr_ko_hfd6', 0.066537),
                    
                    ('igfii_ldlr_ko_nc5', 'igfii_ldlr_ko_nc7', 1.029864),
                    ('igfii_ldlr_ko_hfd5', 'igfii_ldlr_ko_hfd7', 1.243094),
                    
                    ('wt_nc1', 'ldlr_ko_nc4', 1.670599),
                    ('wt_hfd1', 'ldlr_ko_hfd4', 74.018798),
                    ('wt_nc3', 'ldlr_ko_nc6', 0.841111),
                    ('wt_hfd3', 'ldlr_ko_hfd6', 2.805988),
                    
                    ('wt_nc1', 'igfii_ldlr_ko_nc5', 1.314244),
                    ('wt_hfd1', 'igfii_ldlr_ko_hfd5', 1.181026),
                    ('wt_nc3', 'igfii_ldlr_ko_nc7', 0.706487),
                    ('wt_hfd3', 'igfii_ldlr_ko_hfd7', 0.890705),
                    
                    ('ldlr_ko_nc4', 'igfii_ldlr_ko_nc5', 0.773865),
                    ('ldlr_ko_hfd4', 'igfii_ldlr_ko_hfd5', 0.017159),
                    ('ldlr_ko_nc6', 'igfii_ldlr_ko_nc7', 0.843301),
                    ('ldlr_ko_hfd6', 'igfii_ldlr_ko_hfd7', 0.323414),
                   
                    ('ldlr_ko_nc4', 'ldlr_ko_hfd6', 2.850102),
                    ('ldlr_ko_hfd6', 'igfii_ldlr_ko_hfd5', 0.258872),
                    ('wt_hfd1', 'ldlr_ko_hfd6', 4.666616),
                    
                    ('wt_nc1', 'wt_hfd3', 1.805155),
                    ('wt_nc3', 'wt_hfd1', 0.574894),
                    ('wt_nc1', 'wt_hfd2', 1.185841),
                    ('wt_nc3', 'wt_hfd2', 0.609282),
                    ('wt_hfd1', 'wt_hfd2', 1.083545),
                    ('wt_hfd2', 'wt_hfd3', 1.510782),
                   )

# Sets paired by comparison, paired by replicate comparison
# Note that we don't use wt_nc2, wt_hfd2, or ldlr_ko_hfd4 here due to low tag counts.
replicate_sets = ((('wt_nc1', 'wt_hfd1'), ('wt_nc3', 'wt_hfd3')),
                    (('ldlr_ko_nc4', 'ldlr_ko_hfd6'), ('ldlr_ko_nc6', 'ldlr_ko_hfd6')),
                    (('igfii_ldlr_ko_nc5', 'igfii_ldlr_ko_hfd5'), ('igfii_ldlr_ko_nc7', 'igfii_ldlr_ko_hfd7')),
                    
                    (('wt_nc1', 'ldlr_ko_nc4'), ('wt_nc3', 'ldlr_ko_nc6')),
                    (('wt_hfd1', 'ldlr_ko_hfd6'), ('wt_hfd3', 'ldlr_ko_hfd6')),
                    
                    (('wt_nc1', 'igfii_ldlr_ko_nc5'), ('wt_nc3', 'igfii_ldlr_ko_nc7')),
                    (('wt_hfd1', 'igfii_ldlr_ko_hfd5'), ('wt_hfd3', 'igfii_ldlr_ko_hfd7')),
                    
                    (('ldlr_ko_nc4', 'igfii_ldlr_ko_nc5'), ('ldlr_ko_nc6', 'igfii_ldlr_ko_nc7')),
                    (('ldlr_ko_hfd6', 'igfii_ldlr_ko_hfd5'), ('ldlr_ko_hfd6', 'igfii_ldlr_ko_hfd7')),
                   
                   )

pretty_names = {'wt_nc': 'C57Bl6 NC ',
                'wt_hfd': 'C57Bl6 HFD ',
                'ldlr_ko_nc': 'LDLRko NC ',
                'ldlr_ko_hfd': 'LDLRko HFD ',
                'igfii_ldlr_ko_nc': 'IGFII-LDLRko NC ',
                'igfii_ldlr_ko_hfd': 'IGFII-LDLRko HFD ',
                }