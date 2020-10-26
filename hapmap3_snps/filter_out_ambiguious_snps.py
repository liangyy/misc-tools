import numpy as np

SNP_PAIR = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

def filter_out_ambiguious_snps(df):
    '''
    Input a pandas DataFrame with BIM entries: 
    ['chr', 'rsid', 'placeholder', 'pos', 'a1', 'a2']
    '''
    is_ambi = []
    for a, b in zip(list(df.a1), list(df.a2)):
        if SNP_PAIR[a] == b:
            is_ambi.append(True)
        else:
            is_ambi.append(False)
    df = df[np.logical_not(np.array(is_ambi))].reset_index(drop=True)
    return df        

def clean_up_rsid(mylist):
    o = []
    for i in mylist:
        if i.startswith('rs'):
            o.append(i)
        elif i.startswith('SNP'):
            pass
        elif i.startswith('AFFX-SNP'):
            i = i.split('__')[-1]
            o.append(i)
    return o
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='filter_out_ambiguious_snps.py', description='''
        Filter out ambiguious SNPs.
        Will keep SNPs with legitimate rsID only.
        CAUTION: Need to have path to misc-tools/pyutil in the PYTHONPATH!
    ''')
    parser.add_argument('--input', help='''
        A list of SNPs.
    ''')
    parser.add_argument('--input_bim', help='''
        A plink BIM file.
    ''')
    parser.add_argument('--output', help='''
        Output SNP list.
    ''')
    args = parser.parse_args()
 
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    
    import pandas as pd
    from pyutil import load_list, intersection, write_list
    snp_list = load_list(args.input)
    df_bim = pd.read_csv(args.input_bim, sep='\s+', header=None)
    df_bim.columns = ['chr', 'rsid', 'placeholder', 'pos', 'a1', 'a2']
    
    df_bim = filter_out_ambiguious_snps(df_bim)
    snp_list = intersection(snp_list, list(df_bim.rsid))
    
    # before output
    # I found that there are some SNPs with ID like: AFFX-SNP_XXX__rsYYY or SNP_A-ZZZ
    # I cleaned them up by removing all SNPs with SNP_A-ZZZ as ID and 
    # keep rsYYY from AFFX-SNP_XXX__rsYYY 
    snp_list = clean_up_rsid(snp_list)
    
    write_list(snp_list, args.output)
    
