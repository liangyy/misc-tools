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
    for a, b in zip(df.a1.to_list(), df.a2.to_list()):
        if SNP_PAIR[a] == b:
            is_ambi.append(True)
        else:
            is_ambi.append(False)
    df = df[np.logical_not(np.array(df))].reset_index(drop=True)
    return df        
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='filter_out_ambiguious_snps.py', description='''
        Filter out ambiguious SNPs.
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
    
    from pyutil import load_list, intersection, write_list
    snp_list = load_list(args.input)
    df_bim = pd.read_csv(args.input_bim, sep='\s+', header=None)
    df_bim.columns = ['chr', 'rsid', 'placeholder', 'pos', 'a1', 'a2']
    
    df_bim = filter_out_ambiguious_snps(df_bim)
    snp_list = intersection(snp_list, df_bim.rsid.to_list())
    write_list(snp_list, args.output)
    