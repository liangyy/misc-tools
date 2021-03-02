import re
import pandas as pd
import scipy.stats, scipy.io

CHR_WILDCARD = 'CHRNUM'

def replace_str(pat, to_, from_=CHR_WILDCARD):
    return re.sub(from_, str(to_), pat)

def _rs2int(ll):
    oo = - np.ones(length(ll))
    for i, l in enumerate(ll):
        try:
            oo[i] = int(re.sub('rs', '', l))
        except:
            continue
    return oo

def _load_gwas(gwas_file):
    dd = pd.read_parquet(gwas_file)
    dd = dd[['variant_id', 'pval']]
    dd['chisq'] = scipy.stats.chi2.isf(dd.pval, df=1)
    dd['rs_int'] = _rs2int(list(dd.variant_id))
    return dd[['rs_int', 'chisq']]

def load_gwas(gwas_pattern, logging):
    if CHR_WILDCARD in gwas_pattern:
        raise ValueError('Need gwas_pattern and bim_pattern both have or not have {}.'.format(CHR_WILDCARD))
        out = []
        for i in range(1, 23):
            logging.info('load_gwas: chr = {}'.format(i))
            dd = _load_gwas(replace_str(gwas_pattern, i))
            out.append(dd)
        out = pd.concat(out, axis=0).reset_index(drop=True)
    else:
        logging.info('load_gwas: one file')
        out =  _load_gwas(gwas_pattern)  
    return out
            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='format_gwas_to_mat.py', description='''
        Format GWAS into SLD4M format (MATLAB mat file).
        Take GWAS parquet file generated from tensorqtl run 
        (using a customized script).
    ''')
    parser.add_argument('--parquet', help='''
        GWAS parquet file.
        If wildcard CHRNUM is there, we read from chr_num = 1 to 22.
        In the GWAS file, there should be columns:
        variant_id, phenotype_id, pval, b, b_se, maf
    ''')
    parser.add_argument('--output', help='''
        Output MATLAB mat file for SLD4M.
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
    df_gwas = load_gwas(args.parquet, logging)
    n0 = df_gwas.shape[0]
    df_gwas = df_gwas[ df_gwas.rs_int > 0 ].reset_index(drop=True)
    logging.info('Loaded {} SNPs in total and {} have valid rsID.'.format(n0, df_gwas.shape[0]))
    
    logging.info('Writing to disk.')
    # df_gwas.to_csv(args.output, index=False, sep=' ')
    out_dict = {'chisq': dd.chisq.astype(np.float64), 'sumstat_RSIDs_as_ints': dd.rs_int.astype(np.int32)}
    scipy.io.savemat(args.output, out_dict)
    