'''
Created on Oct 2, 2014

@author: karmel

Output TSV for all of the WT and Foxo1 KO ATAC seq samples from the first 
dataset.
'''
from glasslab.dataanalysis.misc.rodrigo.samples import SAMPLES, sample_name
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer
from glasslab.utils.database import get_engine, dataframe_from_query
if __name__ == '__main__':
    yzer = MotifAnalyzer()

    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/' +\
        'Miscellaneous_Collaborations/Rodrigo_CD8s_2014_09/Promoters'
    dirpath = yzer.get_path(dirpath)

    # Get DB engine
    engine = get_engine(
        uri='ec2-23-20-125-153.compute-1.amazonaws.com',
        password='rodrigo#cd8')

    for i, (cond, seq, breed) in enumerate(SAMPLES):
        others = SAMPLES[:i] + SAMPLES[i + 1:]
        sql = '''select distinct on (p1.id) 
    chr."name" as chr_name, p1."start", p1."end", p1.tag_count,
    p1.*,
    '''
        selects, joins = [], []
        for j, (oth_cond, oth_seq, oth_breed) in enumerate(others):
            counter = j + 2
            selects.append('''p{counter}.tag_count as {}{}_{}_tag_count
    '''.format(oth_breed, oth_cond, oth_seq, counter=counter))
            joins.append('''left outer join 
    chipseq.peak_{}cd8tcell_{}_{} p{counter}
    on p1.chromosome_id = p{counter}.chromosome_id
    and p1.start_end && p{counter}.start_end
            '''.format(oth_breed, oth_seq, oth_cond, counter=counter))

        sql += ',\n'.join(selects)
        sql += ''',
    tcf1.tag_count as tcf1_tag_count
    from chipseq.peak_{}cd8tcell_{}_{} p1
    join genome_reference_mm10.chromosome chr 
    on p1.chromosome_id = chr.id
    '''.format(breed, seq, cond)
        sql += ''.join(joins)
        sql += '''left outer join 
    chipseq.peak_cd8tcell_tcf1 tcf1
    on p1.chromosome_id = tcf1.chromosome_id
    and p1.start_end && tcf1.start_end
    left outer join genome_reference_mm10.sequence_transcription_region reg
    on p1.chromosome_id = reg.chromosome_id
    and p1.start_end && int8range(reg.strand*(reg.transcription_end - 100) + 
            abs(reg.strand - 1)*(reg.transcription_start - 100), 
        reg.strand*(reg.transcription_end + 100) + 
            abs(reg.strand - 1)*(reg.transcription_start + 100))
    where reg.id is NOT NULL
    AND p1.tag_count >= 20;
    '''

        # Set up output dir
        sample_prefix = sample_name(cond, seq, breed)
        sample_dirpath = yzer.get_filename(dirpath, sample_prefix)

        # Get data
        data = dataframe_from_query(sql, engine)

        output_file = yzer.get_filename(
            sample_dirpath, sample_prefix + '_promoters.txt')
        data.to_csv(output_file, sep='\t', header=True, index=False)
