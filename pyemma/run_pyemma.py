import pandas as pd

def subset_grm(grm, grm_indiv, target_indiv):
    set_target_indiv = set(target_indiv)
    isin = np.array([ g in set_target_indiv for g in grm_indiv ])
    grm = grm[:, isin][isin, :]
    grm_indiv = list(np.array(grm_indiv)[isin])
    return grm, grm_indiv

def subset_y(df, indiv):
    df_indiv = pd.DataFrame({'indiv': indiv})
    df = pd.merge(df_indiv, df, on='indiv')
    return df.iloc[:, 1:].values

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_pyemma.py', description='''
        Run a series of variance component estimation given GRM and y.
    ''')
    parser.add_argument('--grm', help='''
        Provide the prefix of GRM files. 
        Will assume [prefix].grm.gz and [prefix].grm.id.
    ''')
    parser.add_argument('--reml', action='store_true', help='''
        If specified, it will use reml.
    ''')
    parser.add_argument('--grm_cache', default=None, help='''
        Optional. If specified, will use the cached GRM EVD if file exists.
        If specified but does not exists, will read grm and calculate EVD 
        from scratch and cache the result.
        It is the prefix. Will append reml or mle accordingly. 
        CAUTION: If using GRM cache, the --y_table should have exactly 
        the same set of individuals. 
    ''')
    parser.add_argument('--y_table', nargs='+', help='''
        The table of y along with the column name of individual ID.
    ''')
    parser.add_argument('--y_list', default=None, help='''
        Optional. If specified, will limit the analysis to 
        the phenotypes in the list
    ''')
    parser.add_argument('--output', help='''
        Output result summary in TSV.GZ format.
    ''')
    parser.add_argument('--evd_min_max_ratio', type=float, default=None, help='''
        The cutoff on eigen-vectors and -values in EVD: value / max(value) > evd_min_max_ratio.
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
    
    import numpy as np
    import gzip
    import pickle
    from tqdm import tqdm
    import scipy.stats

    from pyutil import file_exists, read_table, load_list, intersection
    
    import pyemma
    
    if args.grm_cache is not None:
        if args.reml is True:
            grm_cache = args.grm_cache + '.reml.pkl.gz'
        else:
            grm_cache = args.grm_cache + '.mle.pkl.gz'
    else:
        grm_cache = None
    
    logging.info('Loading phenotype table.')
    df_y = read_table(
        args.y_table[0], 
        indiv_col=args.y_table[1]
    )
    pheno_indiv = df_y.indiv.to_list()
    if args.y_list is not None:
        y_list = load_list(args.y_list)
        df_y = df_y[['indiv'] + y_list]
    pheno_list = df_y.columns[1:].to_list()
    
    if grm_cache is None or not file_exists(grm_cache):
        logging.info('Loading GRM from scratch.')
        grm, grm_indiv = pyemma.load_grm(args.grm + '.grm.gz', args.grm + '.grm.id')
        common_indiv = intersection(grm_indiv, pheno_indiv)
        grm, indiv = subset_grm(grm, grm_indiv, common_indiv)
        ymat = subset_y(df_y, indiv)
        
        logging.info('Starting GRM EVD.')
        
        eig_val, eig_vec = pyemma.pyemma_mle_mat_fac(grm, min_max_ratio=args.evd_min_max_ratio)
        to_cache = {'vec': eig_vec, 'val': eig_val, 'indiv': indiv}
        if args.reml is True:
            eig_val_inter, eig_vec_inter = pyemma.pyemma_reml_mat_fac(np.ones((grm.shape[0], 1)), grm, min_max_ratio=args.evd_min_max_ratio)
            to_cache['vec_intercept'] = eig_vec_inter
            to_cache['val_intercept'] = eig_val_inter
        if grm_cache is not None:
            logging.info('Caching GRM EVD.')
            with gzip.open(grm_cache, 'wb') as f:
                pickle.dump(
                    to_cache,
                    f,
                    protocol=4
                )
    else:
        logging.info('Loading GRM from cache.')
        with gzip.open(grm_cache, 'rb') as f:
            tmp = pickle.load(f)
            eig_vec = tmp['vec']
            eig_val = tmp['val']
            if args.reml is True:
                eig_vec_inter = tmp['vec_intercept']
                eig_val_inter = tmp['val_intercept']
            indiv = tmp['indiv']
        ymat = subset_y(df_y, indiv)
        if ymat.shape[0] != len(indiv):
            raise ValueError('We need to have all GRM individaul appear in phenotype table to proceed.')
    
    res = []
    x = np.ones((ymat.shape[0], 1))
    for i in tqdm(range(ymat.shape[1])):
        if args.reml is False:
            res_i = pyemma.pyemma_w_X(ymat[:, i], x, eig_vec, eig_val)
        else:
            res_i = pyemma.pyemma_reml(ymat[:, i], eig_vec_inter, eig_val_inter, eig_vec, eig_val, x)
        res_i['phenotype'] = pheno_list[i]
        res.append(res_i)
    res = pd.concat(res, axis=0)
    
    # use the idea: under the null LR * 2 ~ 0.5 delta(0) + 0.5 * chisq(df=1)
    res['Chisq_pval'] = scipy.stats.chi2.sf(res['LR'] * 2, df=1) * 0.5
    
    res.to_csv(args.output, compression='gzip', sep='\t', index=False)

        
    
    
