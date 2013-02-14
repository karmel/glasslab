'''
Created on Sep 27, 2010

@author: karmel
'''
from django.db import models, connection
from glasslab.genomereference.datatypes import Chromosome, \
    SequencingRun, PeakType
from glasslab.utils.datatypes.basic_model import BoxField
from glasslab.utils.database import execute_query
from glasslab.config import current_settings
from glasslab.sequencing.datatypes.tag import GlassSequencingOutput
   
       
class GlassPeak(GlassSequencingOutput):
    '''
    From MACS::
        
        chr     start   end     length  summit  tags    -10*log10(pvalue)       fold_enrichment
        
    From SICER::
    
        chrom, start, end, ChIP_island_read_count, CONTROL_island_read_count, p_value, fold_change, FDR_threshold
    
    '''
    
    chromosome      = models.ForeignKey(Chromosome)
    start           = models.IntegerField(max_length=12)
    end             = models.IntegerField(max_length=12)
    
    start_end       = BoxField(max_length=255, help_text='This is a placeholder for the PostgreSQL box type.') 
    
    length          = models.IntegerField(max_length=12)
    summit          = models.IntegerField(max_length=12)
    tag_count       = models.DecimalField(max_digits=8, decimal_places=2, null=True, default=None)
    raw_tag_count   = models.DecimalField(max_digits=8, decimal_places=2, null=True, default=None)
    
    score               = models.DecimalField(max_digits=8, decimal_places=2, null=True, default=None)
    p_value             = models.DecimalField(max_digits=6, decimal_places=4, null=True, default=None)
    p_value_exp         = models.IntegerField(max_length=12, null=True, default=None, help_text='Exponent of 10 in p_value x 10^y')
    log_ten_p_value     = models.DecimalField(max_digits=10, decimal_places=2, null=True, default=None)
    fold_enrichment     = models.DecimalField(max_digits=10, decimal_places=2, null=True, default=None)
    fdr_threshold       = models.DecimalField(max_digits=6, decimal_places=4, null=True, default=None)
    fdr_threshold_exp   = models.IntegerField(max_length=12, null=True, default=None, help_text='Exponent of 10 in fdr_threshold x 10^y')
    
    
    @classmethod        
    def create_table(cls, name):
        '''
        Create table that will be used for these peaks,
        dynamically named.
        '''
        cls.set_table_name('peak_' + name)
        
        table_sql = """
        CREATE TABLE "%s" (
            id int4,
            chromosome_id int4,
            "start" int8,
            "end" int8,
            start_end box,
            "length" int4,
            summit int8,
            tag_count decimal(8,2) default NULL,
            raw_tag_count decimal(8,2) default NULL,
            score decimal(8,2) default NULL,
            p_value decimal(6,4) default NULL,
            p_value_exp int4 default NULL,
            log_ten_p_value decimal(10,2) default NULL,
            fold_enrichment decimal(10,2) default NULL,
            fdr_threshold decimal(6,4) default NULL,
            fdr_threshold_exp int4 default NULL
            );
        CREATE SEQUENCE "%s_id_seq"
            START WITH 1
            INCREMENT BY 1
            NO MINVALUE
            NO MAXVALUE
            CACHE 1;
        ALTER SEQUENCE "%s_id_seq" OWNED BY "%s".id;
        ALTER TABLE "%s" ALTER COLUMN id SET DEFAULT nextval('"%s_id_seq"'::regclass);
        ALTER TABLE ONLY "%s" ADD CONSTRAINT %s_pkey PRIMARY KEY (id);
        """ % (cls._meta.db_table, 
               cls._meta.db_table,
               cls._meta.db_table, cls._meta.db_table,
               cls._meta.db_table, cls._meta.db_table,
               cls._meta.db_table, cls.name)
        execute_query(table_sql)
        
        cls.table_created = True
    
    @classmethod
    def add_indices(cls):
        update_query = """
        CREATE INDEX %s_chr_idx ON "%s" USING btree (chromosome_id);
        CREATE INDEX %s_start_end_idx ON "%s" USING gist (start_end);
        """ % (cls.name, cls._meta.db_table,
               cls.name, cls._meta.db_table)
        execute_query(update_query)
            
    @classmethod
    def init_from_macs_row(cls, row):
        '''
        From a standard tab-delimited MACS peak file, create model instance.
        '''
        return cls(chromosome=Chromosome.objects.get(name=str(row[0]).strip()),
                     start=int(row[1]),
                     end=int(row[2]),
                     start_end=(int(row[1]), 0, int(row[2]), 0),
                     length=int(row[3]),
                     summit=int(row[4]),
                     tag_count=str(row[5]),
                     log_ten_p_value=str(row[6]),
                     fold_enrichment=str(row[7])
                     )
    @classmethod
    def init_from_homer_row(cls, row):
        '''
        From a standard tab-delimited Homer peak file, create model instance.
        '''
        connection.close()
        try: p_val = str(row[11]).lower().split('e')
        except IndexError: p_val = None
        return cls(chromosome=Chromosome.objects.get(name=str(row[1]).strip()),
                     start=int(row[2]),
                     end=int(row[3]),
                     start_end=(int(row[2]), 0, int(row[3]), 0),
                     length=int(row[3]) - int(row[2]),
                     tag_count=str(row[5]),
                     score=str(row[7]),
                     p_value=p_val and str(p_val[0]) or None,
                     p_value_exp=p_val and len(p_val) > 1 and p_val[1] or 0,
                     fold_enrichment=str(row[10])
                     )
        
    @classmethod
    def init_from_homer_row_old(cls, row):
        '''
        From a standard tab-delimited Homer peak file, create model instance.
        
        For use with older versions of Homer files.
        '''
        connection.close()
        p_val = str(row[12]).lower().split('e')
        return cls(chromosome=Chromosome.objects.get(name=str(row[1]).strip()),
                     start=int(row[2]),
                     end=int(row[3]),
                     start_end=(int(row[2]), 0, int(row[3]), 0),
                     length=int(row[3]) - int(row[2]),
                     tag_count=str(row[5]),
                     raw_tag_count=str(row[9]),
                     score=str(row[8]),
                     p_value=str(p_val[0]),
                     p_value_exp=len(p_val) > 1 and p_val[1] or 0,
                     fold_enrichment=str(row[11])
                     )
        
        
    @classmethod
    def init_from_sicer_row(cls, row):
        '''
        From a standard tab-delimited SICER peak file, create model instance.
        '''
        p_val = str(row[5]).split('e') # '1.34e-12' --> 1.34, -12
        fdr = str(row[7]).split('e')
        return cls(chromosome=Chromosome.objects.get(name=str(row[0]).strip()),
                     start=int(row[1]),
                     end=int(row[2]),
                     start_end=(int(row[1]), 0, int(row[2]), 0),
                     length=int(row[2]) - int(row[1]),
                     tag_count=str(row[3]),
                     p_value=p_val[0],
                     p_value_exp=len(p_val) > 1 and p_val[1] or 0,
                     fold_enrichment=str(row[6]),
                     fdr_threshold=fdr[0],
                     fdr_threshold_exp=len(fdr) > 1 and fdr[1] or 0
                     )
    @classmethod
    def init_from_bed_row(cls, row):
        '''
        From a standard tab-delimited BED peak file, create model instance.
        '''
        return cls(chromosome=Chromosome.objects.get(name=str(row[0]).strip()),
                     start=int(row[1]),
                     end=int(row[2]),
                     start_end=(int(row[1]), 0, int(row[2]), 0),
                     length=int(row[2]) - int(row[1]),
                     score=str(row[4]),
                     )
        
    @classmethod
    def score_sicer_peaks(cls):
        '''
        Scale tag count to get a density-based score for peaks.
        
        Total tags uses only those that fall into peaks, assuming that 
        the remaining tags are noisy.
        '''
        update_query = """
        UPDATE "%s" SET
            score = (tag_count::numeric/length/
                (SELECT sum(tag_count) FROM "%s"))*(10^10);
        """ % (cls._meta.db_table,
               cls._meta.db_table)
        execute_query(update_query)
    
    _peak_type = None
    @classmethod 
    def peak_type(cls, name=None):
        '''
        Try to determine peak type from table name.
        '''
        if cls._peak_type: return cls._peak_type
        name = (name or cls.name).lower()
        peak_types = PeakType.objects.iterator()
        for p_type in peak_types: 
            if name.find(p_type.type.strip().lower().replace('.','_')) >= 0:
                cls._peak_type = p_type
                break
        if not cls._peak_type: raise Exception('Peak type not found for {0}'.format(name))
        return cls._peak_type
    
    @classmethod 
    def add_record_of_tags(cls, description='', type='ChIP-Seq', peak_type=None, stats_file=None):
        '''
        Add SequencingRun record with the details of this run.
        
        Should be called only after all tags have been added.
        ''' 
        connection.close()
        total_tags, percent_mapped = cls.get_bowtie_stats(stats_file)
        wt, notx, kla, other_conditions, timepoint = cls.parse_attributes_from_name()
        s, created = SequencingRun.objects.get_or_create(source_table=cls._meta.db_table,
                                        defaults={'name': cls.name, 
                                                  'total_tags': total_tags,
                                                  'description': description,
                                                  'cell_type': current_settings.CELL_TYPE,
                                                  'type': type,
                                                  'peak_type': peak_type or cls.peak_type(),
                                                  'percent_mapped': percent_mapped,
                                                  'wt': wt,
                                                  'notx': notx,
                                                  'kla': kla,
                                                  'other_conditions': other_conditions,
                                                  'timepoint': timepoint, }
                                               )
        if not created: 
            s.total_tags = total_tags
            s.percent_mapped = percent_mapped
            s.wt, s.notx, s.kla, s.other_conditions = wt, notx, kla, other_conditions 
            s.timepoint = timepoint  
            s.save() 
        return s
    
class HomerPeak(GlassSequencingOutput):
    '''
    Peaks from Homer output, for easy table loading and comparison.
    '''
    cluster_id      = models.CharField(max_length=255, null=True, blank=True)
    chromosome      = models.ForeignKey(Chromosome)
    strand          = models.IntegerField(max_length=1)
    start           = models.IntegerField(max_length=12)
    end             = models.IntegerField(max_length=12)
    
    start_end       = BoxField(max_length=255, help_text='This is a placeholder for the PostgreSQL box type.') 
    
    peak_score      = models.IntegerField(max_length=12)
    distance_to_tss = models.IntegerField(max_length=12)
    
    annotation          = models.CharField(max_length=255, null=True, blank=True)
    detailed_annotation = models.CharField(max_length=255, null=True, blank=True)
    nearest_promoter    = models.CharField(max_length=100, null=True, blank=True)
    promoter_id         = models.IntegerField(max_length=12, null=True, blank=True)
    nearest_unigene     = models.CharField(max_length=100, null=True, blank=True)
    nearest_refseq      = models.CharField(max_length=100, null=True, blank=True)
    nearest_ensembl     = models.CharField(max_length=100, null=True, blank=True)
    gene_name           = models.CharField(max_length=100, null=True, blank=True)
    gene_alias          = models.CharField(max_length=255, null=True, blank=True)
    gene_description    = models.CharField(max_length=255, null=True, blank=True)

    @classmethod        
    def create_table(cls, name):
        '''
        Create table that will be used for these peaks,
        dynamically named.
        '''
        cls.set_table_name('homer_peak_' + name)
        
        table_sql = """
        CREATE TABLE "%s" (
            id int4,
            cluster_id varchar(255),
            chromosome_id int4,
            strand int2,
            "start" int8,
            "end" int8,
            start_end box,
            "peak_score" int4,
            "distance_to_tss" int4,
            annotation varchar(255),
            detailed_annotation varchar(255),
            nearest_promoter varchar(100),
            promoter_id int4,
            nearest_unigene varchar(100),
            nearest_refseq varchar(100),
            nearest_ensembl varchar(100),
            gene_name varchar(100),
            gene_alias varchar(255),
            gene_description varchar(255)
            );
        CREATE SEQUENCE "%s_id_seq"
            START WITH 1
            INCREMENT BY 1
            NO MINVALUE
            NO MAXVALUE
            CACHE 1;
        ALTER SEQUENCE "%s_id_seq" OWNED BY "%s".id;
        ALTER TABLE "%s" ALTER COLUMN id SET DEFAULT nextval('"%s_id_seq"'::regclass);
        ALTER TABLE ONLY "%s" ADD CONSTRAINT %s_pkey PRIMARY KEY (id);
        """ % (cls._meta.db_table, 
               cls._meta.db_table,
               cls._meta.db_table, cls._meta.db_table,
               cls._meta.db_table, cls._meta.db_table,
               cls._meta.db_table, cls.name)
        execute_query(table_sql)
        
        cls.table_created = True
    
    @classmethod
    def add_indices(cls):
        update_query = """
        CREATE INDEX %s_chr_idx ON "%s" USING btree (chromosome_id);
        CREATE INDEX %s_strand_idx ON "%s" USING btree (strand);
        CREATE INDEX %s_start_end_idx ON "%s" USING gist (start_end);
        """ % (cls.name, cls._meta.db_table,
               cls.name, cls._meta.db_table,
               cls.name, cls._meta.db_table)
        execute_query(update_query)
        
    @classmethod
    def init_from_homer_row(cls, row):
        '''
        From a standard HOMER .xls file row, create a new peak.
        '''
        try:
            return cls(cluster_id=str(row[0]).strip(),
                   chromosome=Chromosome.objects.get(name=str(row[1]).strip()),
                     strand=str(row[4]).strip() == '-' and 1 or 0,
                     start=row[2] and int(row[2]) or None,
                     end=row[3] and int(row[3]) or None,
                     start_end=(int(row[2]), 0, int(row[3]), 0),
                     peak_score=row[5] and int(row[5]) or None,
                     annotation=str(row[7]).strip(),
                     detailed_annotation=str(row[8]).strip(),
                     distance_to_tss=row[9] and int(row[9]) or None,
                     nearest_promoter=str(row[10]).strip(),
                     promoter_id=row[11] and int(row[11]) or None,
                     nearest_unigene=str(row[12]).strip(),
                     nearest_refseq=str(row[13]).strip(),
                     nearest_ensembl=str(row[14]).strip(),
                     gene_name=str(row[15]).strip(),
                     gene_alias=str(row[16]).strip(),
                     gene_description=str(row[17]).strip()
                     )
        except Exception: 
            print row
            raise