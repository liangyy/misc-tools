if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='gen_lookup_table.py', description='''
        Generate lookup table for HapMap SNP list.
        CAUTION: Need to have path to 
        misc-tools/liftover_snp and misc-tools/pyutil
        in the PYTHONPATH!
    ''')
    parser.add_argument('--input', help='''
        The list of HapMap SNPs.
    ''')
    parser.add_argument('--input_map', help='''
        The HapMap plink MAP file.
    ''')
    parser.add_argument('--liftover_chain', default=None, help='''
        If want to liftover to another genome build, specify the chain file.
    ''')
    parser.add_argument('--output', help='''
        Output SNP lookup table.
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
    import numpy as np
    import lib as lo_lib
    from pyutil import load_list
    
    
    logging.info('Loading SNP list.')
    snp_list = load_list(args.input)
    
    logging.info('Loading SNP MAP file.')
    df_bim = pd.read_csv(args.input_map, sep='\s+', header=None)
    df_bim.columns = ['chr', 'rsid', 'placeholder', 'pos']
    
    logging.info('Intersecting list and MAP.')
    df_bim = df_bim[ df_bim.rsid.isin(snp_list) ].reset_index(drop=True)
    df_bim.chr = df_bim.chr.astype(str)
    # df_bim.pos = df_bim.pos.astype(str)
    logging.info('There are {} SNPs left.'.format(df_bim.shape[0]))
    
    if args.liftover_chain is not None:
        logging.info('Doing liftover, chain file = {}'.format(args.liftover_chain))
        df_lifted = lo_lib.liftover(df_bim.chr, df_bim.pos, args.liftover_chain)
        df_bim.pos = df_lifted.liftover_pos
        df_bim.chr = df_lifted.liftover_chr
        df_bim = df_bim[ ~ df_bim.pos.isna() ].reset_index(drop=True)
        logging.info('There are {} SNPs left after liftover.'.format(df_bim.shape[0]))
    
    df_bim['start'] = df_bim.pos
    df_bim.rename(columns={'pos': 'end', 'rsid': 'name'}, inplace=True)
    df_bim = df_bim[['chr', 'start', 'end', 'name']]
    
    logging.info('Writing to disk.')
    df_bim.to_csv(args.output, sep='\t', compression='gzip', index=False)
    



