'''
Created on Oct 2, 2014

@author: karmel

Output TSV for all of the WT and Foxo1 KO ATAC seq samples from the LCMV
timecourse.

If running on Mac, make sure psycopg2 is set up:
sudo ln -s /Library/PostgreSQL/9.2/lib/libssl.1.0.0.dylib /usr/lib
sudo ln -s /Library/PostgreSQL/9.2/lib/libcrypto.1.0.0.dylib /usr/lib

'''
from glasslab.dataanalysis.motifs.motif_analyzer import MotifAnalyzer
from glasslab.utils.database import get_engine, dataframe_from_query
if __name__ == '__main__':
    yzer = MotifAnalyzer()

    dirpath = 'karmel/Desktop/Projects/GlassLab/Notes_and_Reports/' +\
        'Miscellaneous_Collaborations/Rodrigo_CD8s_2014_09/Enhancers_set2'
    dirpath = yzer.get_path(dirpath)

    # Get DB engine
    engine = get_engine(
        uri='ec2-23-20-125-153.compute-1.amazonaws.com',
        password='rodrigo#cd8')

    # Set up names
    timepoints = (0, 7, 12, 35)
    celltypes = ('klrghi', 'klrglo')
    breeds = ('', 'foxo1_ko_')

    samples = []
    short_names = []
    for breed in breeds:
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

    # Run and save output from sql queries
    for i, sample in enumerate(samples):
        curr_name = short_names[i]
        others = samples[:i] + samples[i + 1:]
        oth_names = short_names[:i] + short_names[i + 1:]
        sql = '''-- {}
    select distinct on (p1.id)
    chr."name" as chr_name, p1."start", p1."end", p1.tag_count,
    p1.*,
    '''.format(sample)
        selects, joins = [], []
        for j, other_sample in enumerate(others):
            counter = j + 2
            selects.append(
                '''p{counter}.tag_count as {}_tag_count'''.format(
                    oth_names[j], counter=counter))
            joins.append('''
    left outer join
    chipseq.peak_{} p{counter}
    on p1.chromosome_id = p{counter}.chromosome_id
    and p1.start_end && p{counter}.start_end
    '''.format(other_sample, counter=counter))

        sql += ',\n'.join(selects)
        sql += '''
    from chipseq.peak_{} p1
    join genome_reference_mm10.chromosome chr 
    on p1.chromosome_id = chr.id
    '''.format(sample)
        sql += ''.join(joins)
        sql += '''
    left outer join genome_reference_mm10.sequence_transcription_region reg
    on p1.chromosome_id = reg.chromosome_id
    and p1.start_end && reg.start_site_1000
    where reg.id is NULL;
    '''
        print(sql)
        # Set up output dir
        sample_path = yzer.get_and_create_path(dirpath, curr_name)

        # Get data
        data = dataframe_from_query(sql, engine)

        output_file = yzer.get_path(sample_path, curr_name + '_enhancers.txt')
        data.to_csv(output_file, sep='\t', header=True, index=False)
