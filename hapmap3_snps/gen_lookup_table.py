def gen_snpid(chr, pos, a1, a2, build):
    res = []
    for x1, x2, x3, x4 in zip(chr, pos, a1, a2):
        res.append(f'chr{x1}_{x2}_{x3}_{x4}_{build}')
    return res

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
    parser.add_argument('--input_bim', help='''
        The HapMap plink BIM file.
    ''')
    parser.add_argument('--input_frq', help='''
        The HapMap plink FRQ file.
    ''')
    parser.add_argument('--liftover_chain', default=None, help='''
        If want to liftover to another genome build, specify the chain file. 
        And target build tag.
    ''')
    parser.add_argument('--target_build', help='''
        The target genome build of output metadata.
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
    import lib as lo_lib
    from pyutil import load_list
    
    
    logging.info('Loading SNP list.')
    snp_list = load_list(args.input)
    
    logging.info('Loading SNP BIM file.')
    df_map = pd.read_csv(args.input_bim, sep='\s+', header=None)
    df_bim.columns = ['chr', 'rsid', 'placeholder', 'pos', 'a1', 'a2']
    
    logging.info('Loading SNP FRQ file.')
    df_frq = pd.read_csv(args.input_frq, sep='\s+', header=None)
    df_bim = pd.merge(df_bim, df_frq[['CHR', 'SNP', 'MAF']], left_on=['chr', 'rsid'], right_on=['CHR', 'SNP'])
    
    logging.info('Intersecting list and MAP.')
    df_bim = df_bim[ df_bim.rsid.isin(snp_list) ].reset_index(drop=True)
    df_bim.chr = df_bim.chr.astype(str)
    # df_bim.pos = df_bim.pos.astype(str)
    logging.info('There are {} SNPs left.'.format(df_bim.shape[0]))
    
    if args.liftover_chain is not None:
        logging.info('Doing liftover, chain file = {}'.format(args.liftover_chain))
        df_lifted = lo_lib.liftover(list(df_bim.chr), list(df_bim.pos), args.liftover_chain)
        df_bim.pos = df_lifted.liftover_pos
        df_bim.chr = df_lifted.liftover_chr
        df_bim = df_bim[ ~ df_bim.pos.isna() ].reset_index(drop=True)
        logging.info('There are {} SNPs left after liftover.'.format(df_bim.shape[0]))
    
    df_bim.rename(
        columns={
            'pos': 'position', 'chr': 'chromosome', 
            'a2': 'allele_0', 'a1': 'allele_1', 
            'allele_1_frequency': 'MAF'
        }, 
        inplace=True
    )
    df_bim['chromosome'] = df_bim.chromosome.apply(lambda x: re.sub('chr', '', str(x)))
    df_bim['id'] = gen_snpid(df_bim.chromosome, df_bim.position, df_bim.allele_0, df_bim.allele_1, args.target_build)
    df_bim = df_bim[['chromosome', 'position', 'id', 'allele_0', 'allele_1', 'allele_1_frequency', 'rsid']]
    
    logging.info('Writing to disk.')
    df_bim.to_csv(args.output, sep='\t', compression='gzip', index=False)
    



