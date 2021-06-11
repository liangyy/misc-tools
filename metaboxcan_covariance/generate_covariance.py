from collections import OrderedDict
import numpy as np

def flip_geno_by_direction(geno, direction):
    vec2 = np.ones((1, geno.shape[1])) * 2
    dire = direction[np.newaxis, :]
    geno = vec2 * (dire == -1) + geno * dire
    return geno

def snp_isin_data(snp, cov_dict):
    if snp.chrom in cov_dict:
        df_snp = cov_dict[snp.chrom][1]
        idx = df_snp[ df_snp.snpid == snp.snpid ].idx.tolist()
        if len(idx) == 0:
            return False, 'NA'
        elif len(idx) == 1:
            return True, idx[0]
        # return snp.snpid in cov_dict[snp.chrom][1].snpid.tolist(),  
    else:
        return False, 'NA'

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='generate_covariance.py', description='''
        Prepare genotype covariance for MetaboXcan PredictDBs.
        We take an all-in-one genotype file in PLINK BED format.
        Genotype IO is transethnic_prs.util.genotype_io.PlinkBedIO.
        PredictDB IO is SPrediXcan2PTRS.util.db.WeightDB.
        WARNING: we don't handle ambiguious SNPs and cases where multiple SNPs 
        have the same genomic position. 
    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--genotype_bed', help='''
        Path to the genotype BED file. 
    ''')
    parser.add_argument('--predictdb', help='''
        PredictDB file. 
    ''')
    parser.add_argument('--snpid_to_report', nargs='+', default=['rsid'], help='''
        rsid or varID.
    ''')
    # parser.add_argument('--mode', nargs='+', help='''
    #     Indicate the mode and parameter of the mode for genotype covariance:
    #     1. banded [band-size]; 
    #     2. cap [threshold-of-cap];
    #     3. naive [f32|f64];
    #     4. evd [min-max-threshold]. 
    # ''')
    parser.add_argument('--output_prefix', help='''
        Output covariance in TXT.GZ format.
        Will add extension depending on --snpid_to_report. 
        (ready to go with SPrediXcan.py BUT DO NOT use `streamline` mode).
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
    
    import gzip
    from tqdm import tqdm
    import pandas as pd
    import SPrediXcan2PTRS.util.db as db
    import transethnic_prs.util.genotype_io as genoio
    
    logging.info('Loading SNP ID to report.')
    id_to_report = []
    for i in args.snpid_to_report:
        if i not in ['rsid', 'varID']:
            raise ValueError(f'args.snpid_to_report can only take rsid or varID but {i} is given.')
        id_to_report.append(i)
    
    logging.info('Loading predictdb.')
    func_get_snp_meta = lambda x: db.get_snp_meta_from_varID(x)
    weight = db.WeightDB(args.predictdb, func_get_snp_meta)

    logging.info('Loading genotype.')
    geno = genoio.PlinkBedIO(args.genotype_bed)
    
    out_dict = OrderedDict()
    for i in id_to_report:
        fn = '{}.{}.txt.gz'.format(args.output_prefix, i)
        f = gzip.open(fn, 'wt')
        f.write('GENE RSID1 RSID2 VALUE\n')
        out_dict[i] = f

    genes = weight.get_gene_info()
    logging.info('Going through one gene at a time: ngene = {}'.format(len(genes)))
    for gidx in tqdm(range(len(genes))):
        gene = genes[gidx]
        snps = weight.get_variants_of_gene([gene])
        # annotate with snp ids in weight file
        snps = pd.merge(
            snps, 
            weight.df_var[['chrom', 'position', 'effect_allele', 'non_effect_allele', 'rsid', 'varID']], 
            how='left', 
            on=['chrom', 'position', 'effect_allele', 'non_effect_allele']
        )
        ndup = snps[['chrom', 'position']].duplicated().sum()
        if ndup > 0:
            raise ValueError(f'Gene = {gene} has multiple SNPs sharing the same genomic position. Exit!')
        # annotate with snpids
        snps['idx'] = [ i for i in range(snps.shape[0]) ]
        snpids = genoio.snpinfo_to_snpid(snps.chrom, snps.position, snps.effect_allele, snps.non_effect_allele)
        snps = pd.merge(snps, snpids[['idx', 'direction', 'snpid']], on='idx', how='left')
        snps.drop(columns=['idx'], inplace=True)
        cov_dict = OrderedDict()
        for cc in snps.chrom.unique():
            snps_cc = snps[ snps.chrom == cc ].reset_index(drop=True).copy()
            if snps_cc.snpid.isna().sum() > 0:
                logging.info(f'Gene = {gene}: It is likely that there are ambiguious SNPs. They will be discarded.')
            snpid = snps_cc.snpid[~snps_cc.snpid.isna()].values
            if len(snpid) == 0:
                continue
            mat, snplist = geno.load(snpid, return_snplist=True)
            if mat is None:
                continue
            df_snp = pd.merge(snplist[['snpid']], snps_cc[['snpid', 'rsid', 'varID', 'direction', 'chrom']], how='left', on='snpid')
            df_snp['idx'] = [ i for i in range(df_snp.shape[0]) ]
            mat = flip_geno_by_direction(mat, df_snp.direction.values)
            cov_ = np.cov(mat.T)
            if len(cov_.shape) == 0:
                cov_ = cov_[np.newaxis, np.newaxis]
            cov_dict[cc] = (cov_, df_snp)
        # write to file
        for i in range(snps.shape[0]):
            snpi = snps.iloc[i]
            snpi_isin_data, idxi = snp_isin_data(snpi, cov_dict)  
            for j in range(i + 1):
                snpj = snps.iloc[j]
                snpj_isin_data, idxj = snp_isin_data(snpj, cov_dict)  # snpj.snpid in cov_dict[snpj.chrom][1].snpid
                if snpi_isin_data is False or snpj_isin_data is False:
                    val = 'NA'
                else:
                    if snpi.chrom == snpj.chrom:
                        cov_this = cov_dict[snpi.chrom][0]
                        val = cov_this[idxi, idxj]
                    else:
                        val = 0
                for k, v in out_dict.items():
                    v.write('{} {} {} {}\n'.format(
                        gene, snpi[k], snpj[k], val
                    ))
    for v in out_dict.values():
        v.close()
    
    
    logging.info('Done.') 
