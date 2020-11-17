import pandas as pd

def transpose_df(df, col):
    return df.set_index(col).T

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_trans_qtl.py', description='''
        Run GWAS of quantitative traits. 
        Need to export PYTHONPATH=path-to/misc-tools/pyutil:path-to/tensorqtl
    ''')
    parser.add_argument('--geno_bed_prefix', help='''
        Plink binary PED format (plink1).
    ''')
    parser.add_argument('--phenotype_table', nargs='+', help='''
        In parquet or csv format.
        Specify the filename followed by the column of individual ID.
    ''')
    parser.add_argument('--covariate_table', nargs='+', help='''
        In parquet or csv format.
        Specify the filename followed by the column of individual ID.
    ''')
    parser.add_argument('--covariate_yaml', default=None, help='''
        Specify the list of covariate to use and 
        the type of the covariate (continuous or categorical).
    ''')
    parser.add_argument('--output_prefix', help='''
        Will output parquet file for each phenotype.
        [output_prefix].[phenotype_name].parquet
    ''')
    parser.add_argument('--map_trans_params', default=None, help='''
        A YAML file containing the arguments to feed into tensorqtl.map_trans.
    ''')
    parser.add_argument('--individual_list', help='''
        The list of individuals to include in the analysis.
    ''')
    parser.add_argument('--individual_list_exclude', default=None, help='''
        The list of individuals to exclude from the analysis.
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
    from tqdm import tqdm
    from pyutil import load_list
    from util import exclude_b_from_a
    import genotypeio
    from tensorqtl import map_trans
    
    logging.info('Loading tables.')
    df_covar = load_covariate(args.covariate_table, args.covariate_yaml)
    df_pheno = load_phenotype(args.phenotype_table)
    indiv_lists = [df_covar.indiv.to_list(), df_pheno.indiv.to_list()]
    if args.individual_list is not None:
        indiv_lists.append(load_list(args.individual_list))
    indiv_list = take_intersect(indiv_lists)
    
    if args.individual_list_exclude is not None:
        indiv_list = exclude_b_from_a(
            a=indiv_list, 
            b=load_list(args.individual_list_exclude)
        )
    
    indiv_list = sorted(indiv_list)
    
    df_covar = rearrange_rows(df_covar, indiv_list)
    df_pheno = rearrange_rows(df_pheno, indiv_list)
    df_pheno = transpose_df(df_pheno, col='indiv')
    df_covar.set_index('indiv', inplace=True)
    logging.info('There are {} individauls being included.'.format(df_covar.shape[0]))
    
    logging.info('Loading genotypes.')
    pr = genotypeio.PlinkReader(
        args.geno_bed_prefix, 
        select_samples=df_pheno.columns, 
        dtype=np.int8
    )
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    genotype_df = pd.DataFrame(pr.load_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
    
    if args.map_trans_params is None:
        map_args = {}
    else:
        map_args = read_yaml(args.map_trans_params)
    logging.info('Arguments used for tensorqtl: ', map_args)
    
    logging.info('Start mapping.')
    pairs_df = trans.map_trans(
        genotype_df, df_pheno, df_covar, 
        logger=None, **map_args
    )
    
    logging.info('Finish mapping and writing to disk.')
    pheno_list = list(df_pheno.index)
    for pheno in tqdm(pheno_list):
        sub = pairs_df[pairs_df.phenotype_id == pheno].reset_index(drop=True)
        if sub.shape[0] > 0:
            sub.to_parquet('{}.{}.parquet'.format(args.output_prefix, pheno), index=False)
    logging.info('Done.')
    