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
            following columns: variant_id, chromosome, position (1-based). 
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
    It has the use_cols columns listing chromosome, start, end (1-based) in order.
    It will intersect with the bedfile and 
    annotate with the columns listed in indexes_in_bed (1-based). 
    Return all rows in df_region regardless of being intersected or not.
    '''
    tmp_prefix = clean_prefix(tmp_prefix)
    
    if len(use_cols) != 3:
        raise ValueError('Need three columns for use_cols.')
        
    df_region2bed = df_region[use_cols].copy()
    df_region2bed.columns = ['chromosome', 'start', 'end']
    df_region2bed.start = df_region2beds.start - 1  # change to 0-based
    df_region2bed['region_identifier'] = [ i for i in range(df_region2bed.shape[0]) ]
    save_bed(df_region2bed, f'{tmp_prefix}.bed.gz')
    
    sys_call = f'bedtools intersect -a {tmp_prefix}.bed.gz -b {annot_bed} -wa -wb | gzip > {tmp_prefix}.join.bed.gz'
    os.system(sys_call)
    ee = pd.read_csv(f'{tmp_prefix}.join.bed.gz', compression='gzip', sep='\t', header=None)
    
    cols_to_keep = [ i for i in range(df_region2bed.shape[1]) ] + [ df_region2bed.shape[1] + i - 1 for i in indexes_in_bed ]
    names_to_use = ['chromosome', 'start', 'end', 'region_identifier'] + [ 'annot_{}' for i in range(len(indexes_in_bed)) ]
    ee = ee.iloc[:, cols_to_keep]
    ee.columns = names_to_use
    
    df_return = df_region.copy()
    df_return = pd.merge(df_return, ee.iloc[:, 3:], on='region_identifier', how='left')
    df_return = df_return.drop(columns='region_identifier')
    
    return df_return 
    
def annotate_region_with_df(df_region, df2, use_cols, df2_region_cols, df2_other_cols, suffix, tmp_prefix):
    '''
    Annotate df with df2.
    df_region: pandas DataFrame containing information of region in each row. 
    It has the use_cols columns listing chromosome, start, end (1-based) in order.
    df2 should have chromosome, start, end (1-based) in order 
    which are listed in df2_region_cols.
    Specify the columns to add in df2_other_cols 
    '''
    tmp_prefix = clean_prefix(tmp_prefix)
    
    # save df2 as bed file
    if len(df2_region_cols) != 3:
        raise ValueError('Need three columns for use_cols.')
    df2_to_bed = df2[df2_region_cols + df2_other_cols].copy()
    df2_to_bed[df2_region_cols[1]] = df2_to_bed[df2_region_cols[1]] - 1  # change to 0-based
    bed_df2 = f'{tmp_prefix}.df2.bed.gz'
    save_bed(df2_to_bed, bed_df2)
    
    indexes_in_bed = [ 3 + i for i in range(len(df2_other_cols)) ]
    tmp = annotate_region_with_bed(df_region, bed_df2, use_cols, indexes_in_bed)
    tmp.columns[-len(df2_other_cols):] = [ f'{i}_{suffix}' for i in df2_other_cols ]
    
    return tmp
    