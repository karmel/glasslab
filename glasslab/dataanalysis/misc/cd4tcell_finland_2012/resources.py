'''
Created on Nov 7, 2012

@author: karmel
'''
comparison_sets = (('wt_nc1', 'wt_nc3', 2.056744),
('ldlr_ko_hfd6', 'igfii_ldlr_ko_hfd7', 1.708535),
('wt_nc1', 'ldlr_ko_nc4', 1.671835),
('wt_hfd1', 'wt_hfd3', 1.666612),
('wt_hfd2', 'wt_hfd3', 1.512430),
('ldlr_ko_hfd6', 'igfii_ldlr_ko_hfd5', 1.361794),
('wt_nc1', 'igfii_ldlr_ko_nc5', 1.332319),
('wt_nc1', 'igfii_ldlr_ko_hfd5', 1.282950),
('igfii_ldlr_ko_hfd5', 'igfii_ldlr_ko_hfd7', 1.237479),
('igfii_ldlr_ko_nc7', 'igfii_ldlr_ko_hfd7', 1.196474),
('wt_nc1', 'wt_hfd2', 1.185470),
('wt_hfd1', 'igfii_ldlr_ko_hfd5', 1.181884),
('wt_nc1', 'wt_hfd1', 1.088081),
('wt_hfd1', 'wt_hfd2', 1.082038),
('igfii_ldlr_ko_nc5', 'igfii_ldlr_ko_nc7', 1.027743),
('igfii_ldlr_ko_nc5', 'igfii_ldlr_ko_hfd5', 0.982320),
('ldlr_ko_nc4', 'ldlr_ko_nc6', 0.948171),
('wt_nc1', 'ldlr_ko_hfd6', 0.926051),
('wt_nc3', 'wt_hfd3', 0.923376),
('wt_hfd3', 'igfii_ldlr_ko_hfd7', 0.889828),
('wt_hfd1', 'ldlr_ko_hfd6', 0.854582),
('ldlr_ko_nc6', 'igfii_ldlr_ko_nc7', 0.845681),
('wt_nc3', 'ldlr_ko_nc6', 0.842199),
('wt_nc3', 'igfii_ldlr_ko_hfd7', 0.830561),
('ldlr_ko_nc4', 'igfii_ldlr_ko_nc5', 0.775391),
('wt_nc3', 'igfii_ldlr_ko_nc7', 0.701587),
('wt_nc3', 'wt_hfd2', 0.605478),
('ldlr_ko_nc6', 'ldlr_ko_hfd6', 0.579883),
('wt_nc3', 'wt_hfd1', 0.574210),
('ldlr_ko_nc4', 'ldlr_ko_hfd6', 0.549552),
('wt_hfd3', 'ldlr_ko_hfd6', 0.507957),
('wt_nc3', 'ldlr_ko_hfd6', 0.480735),
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
                   
                    (('wt_nc1', 'ldlr_ko_hfd6'), ('wt_nc3', 'ldlr_ko_hfd6')),
                    (('wt_nc1', 'igfii_ldlr_ko_hfd5'), ('wt_nc3', 'igfii_ldlr_ko_hfd7')),
                    
                   )

pretty_names = {'wt_nc': 'C57Bl6 NC ',
                'wt_hfd': 'C57Bl6 HFD ',
                'ldlr_ko_nc': 'LDLRko NC ',
                'ldlr_ko_hfd': 'LDLRko HFD ',
                'igfii_ldlr_ko_nc': 'IGFII-LDLRko NC ',
                'igfii_ldlr_ko_hfd': 'IGFII-LDLRko HFD ',
                }