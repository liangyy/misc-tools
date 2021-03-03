import re
import pandas as pd
import numpy as np
import scipy.stats, scipy.io

CHR_WILDCARD = 'CHRNUM'

def replace_str(pat, to_, from_=CHR_WILDCARD):
    return re.sub(from_, str(to_), pat)

def _rs2int(ll):
    oo = - np.ones(len(ll))
    for i, l in enumerate(ll):
        try:
            oo[i] = int(re.sub('rs', '', l))
        except:
            continue
    return oo

def _to_dict(ll):
    odict = {}
    for l in ll:
        v, k = l.split(':')
        odict[k] = v
    return odict

def _load_gwas(gwas_file, extra=None):
    if extra is None:
        dd = pd.read_parquet(gwas_file)
    else:
        dd = pd.read_csv(gwas_file, compression='gzip', sep='\s+')
        dd.rename(columns=_to_dict(extra), inplace=True)
    dd['chisq'] = (dd.b / dd.b_se) ** 2
    dd['rs_int'] = _rs2int(list(dd.variant_id))
    return dd[['rs_int', 'chisq']]

def load_gwas(gwas_pattern, logging, extra=None):
    if CHR_WILDCARD in gwas_pattern:
        out = []
        for i in range(1, 23):
            logging.info('load_gwas: chr = {}'.format(i))
            dd = _load_gwas(replace_str(gwas_pattern, i), extra=extra)
            out.append(dd)
        out = pd.concat(out, axis=0).reset_index(drop=True)
    else:
        logging.info('load_gwas: one file')
        out =  _load_gwas(gwas_pattern, extra=extra)  
    return out
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='format_gwas_to_mat.py', description='''
        Format GWAS into SLD4M format (MATLAB mat file).
        Take GWAS parquet file generated from tensorqtl run 
        (using a customized script).
    ''')
    parser.add_argument('--gwas_file', help='''
        GWAS parquet file.
        If wildcard CHRNUM is there, we read from chr_num = 1 to 22.
        In the GWAS file, there should be columns:
        variant_id, phenotype_id, pval, b, b_se, maf
    ''')
    parser.add_argument('--output', help='''
        Output MATLAB mat file for SLD4M.
    ''')
    parser.add_argument('--as_text_gz', default=None, nargs='+', help='''
        If GWAS file is text gz. Set it and specify columns:
        variant_id:SNP b:Beta b_se:se
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
    df_gwas = load_gwas(args.gwas_file, logging, extra=args.as_text_gz)
    n0 = df_gwas.shape[0]
    df_gwas = df_gwas[ df_gwas.rs_int > 0 ].reset_index(drop=True)
    logging.info('Loaded {} SNPs in total and {} have valid rsID.'.format(n0, df_gwas.shape[0]))
    
    logging.info('Writing to disk.')
    # df_gwas.to_csv(args.output, index=False, sep=' ')
    out_dict = {'chisq': df_gwas.chisq.astype(np.float64).values[:, np.newaxis], 'sumstat_RSIDs_as_ints': df_gwas.rs_int.astype(np.int32).values[:, np.newaxis]}
    scipy.io.savemat(args.output, out_dict)
    
