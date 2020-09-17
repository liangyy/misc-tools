import os
import pandas as pd
import ntpath

def clean_prefix(ff):
    '''
    Make sure that the prefix ff is not longer than 240 characters
    '''
    basename = ntpath.basename(ff)
    if len(basename) > 240:
        new_name = basename[:240]
        dir_ = os.path.dirname(ff)
        return f'{dir_}/{new_name}'
    else:
        return ff

def save_bed(ff, filename):
    ff = ff.sort_values(by=['chromosome', 'start'], ascending=True, na_position='first')
    ff.to_csv(filename, sep='\t', index=False, header=False, compression='gzip')

def snp2bed(df_snp):
    '''
    df_snp: pandas DataFrame containing SNP information in each row. It has the
            following columns: variant_id, chromosome, position.
    Return: a pandas DataFrame with chromosome, start, end, variant_id (in BED format)
    '''
    df_snp2bed = df_snp.copy()
    df_snp2bed['start'] = df_snp['position']
    df_snp2bed['end'] = df_snp['position'] + 1
    df_snp2bed = df_snp2bed[['chromosome', 'start', 'end', 'variant_id']]
    return df_snp2bed

def intersect_with_bed(df_snp, annot_bed, inplace=True, tmp_prefix='test'):
    '''
    df_snp: pandas DataFrame containing SNP information in each row. It has the
            following columns: variant_id, chromosome, position.
    annot_bed: the bed file to work with.
    tmp_prefix: the function will generate some intermediate files using this prefix. 
                Make sure that there is no more than one function call using this 
                prefix at the same time.
    '''
    tmp_prefix = clean_prefix(tmp_prefix)
    df_snp2bed = snp2bed(df_snp)
    save_bed(df_snp2bed, f'{tmp_prefix}.bed.gz')
    sys_call = f'bedtools intersect -a {tmp_prefix}.bed.gz -b {annot_bed} | gzip > {tmp_prefix}.join.bed.gz'
    os.system(sys_call)
    ee = pd.read_csv(f'{tmp_prefix}.join.bed.gz', compression='gzip', sep='\t', header=None)
    annot = df_snp.variant_id.isin(ee.iloc[:, 3]) * 1
    os.remove(f'{tmp_prefix}.join.bed.gz')
    os.remove(f'{tmp_prefix}.bed.gz')
    if inplace is True:
        df_snp['annot'] = annot
    else:
        df_snp2 = df_snp.copy()
        df_snp2['annot'] = annot
        return df_snp2
    
def annotate_region_with_bed(df_region, bedfile, use_cols, indexes_in_bed, tmp_prefix='test'):
    '''
    df_region: pandas DataFrame containing information of region in each row. 
    It has the use_cols columns listing chromosome, start, end in order.
    It will intersect with the bedfile and 
    annotate with the columns listed in indexes_in_bed. 
    '''
