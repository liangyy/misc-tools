import pandas as pd
import re

CHR_WILDCARD = 'CHRNUM'

def replace_str(pat, to_, from_=CHR_WILDCARD):
    return re.sub(from_, str(to_), pat)

def read_bim(fn):
    tmp = pd.read_csv(fn, sep='\s+', header=None)
    tmp.columns = ['chr', 'SNP', 'PLACEHOLDER', 'pos', 'A2', 'A1']
    return tmp

def _load_gwas(gwas_file, bim_file):
    dd = pd.read_parquet(gwas_file)
    bb = read_bim(bim_file)
    dd = pd.merge(dd, bb, left_on='variant_id', right_on='SNP')
    dd = dd[['SNP', 'A1', 'A2', 'maf', 'b', 'b_se', 'pval']]
    dd.rename(
        columns={'maf': 'freq', 'b_se': 'se', 'pval': 'p'},
        inplace=True
    )
    return dd

def load_gwas(gwas_pattern, bim_pattern, logging):
    if CHR_WILDCARD in gwas_pattern:
        if '' not in bim_pattern:
            raise ValueError('Need gwas_pattern and bim_pattern both have or not have {}.'.format(CHR_WILDCARD))
        out = []
        for i in range(1, 23):
            logging.info('load_gwas: chr = {}'.format(i))
            dd = _load_gwas(
                replace_str(gwas_pattern, i),
                replace_str(bim_pattern, i)
            )
            out.append(dd)
        out = pd.concat(out, axis=0).reset_index(drop=True)
    else:
        logging.info('load_gwas: one file')
        out =  _load_gwas(gwas_pattern, bim_pattern)  
    return out
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='format_gwas.py', description='''
        Format GWAS into GCTB format (MA file).
        Take GWAS parquet file generated from tensorqtl run 
        (using a customized script).
    ''')
    parser.add_argument('--parquet', help='''
        GWAS parquet file.
        If wildcard {chr_num} is there, we read from chr_num = 1 to 22.
        In the GWAS file, there should be columns:
        variant_id, phenotype_id, pval, b, b_se, maf
    ''')
    parser.add_argument('--geno_bim', help='''
        PLINK BIM file to extract the allele information 
        (use the same BIM file as the one used for tensorqtl).
    ''')
    parser.add_argument('--gwas_sample_size', type=int, help='''
        GWAS sample size.
    ''')
    parser.add_argument('--output', help='''
        Output GCTB MA file.
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

    logging.info('Loading GWAS.')
    df_gwas = load_gwas(args.parquet, args.geno_bim, logging)
    df_gwas['n'] = args.gwas_sample_size
    
    logging.info('Writing to disk.')
    df_gwas.to_csv(args.output, index=False, sep=' ')
