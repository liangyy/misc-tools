import sqlite3
import pandas as pd

PAIRED_BASE = {
    'A': 'T', 'T': 'A',
    'G': 'C', 'C': 'G'
}

def read_list(fn):
    o = []
    with open(fn, 'r') as f:
        for i in f:
            o.append(i.strip())
    return o

def write_list(fn, mylist):
    with open(fn, 'w') as f:
        for row in mylist:
            f.write(row + '\n')

def extract_var_df(bgi, snpset):
    '''
    Extract from BGI the list of SNPs with rsID in SNP set (snpset).
    Drop the rows with duplicated rsIDs.
    Annotate with whether it is ambiguous SNPs.
    '''
    with sqlite3.connect(bgi) as conn:
        variants = conn.execute('select * from Variant').fetchall()
    o = {'rsid': [], 'a0': [], 'a1': []}
    for _, _, rsid, _, a0, a1, _, _ in variants:
        if rsid in snpset:
            o['rsid'].append(rsid)
            o['a0'].append(a0)
            o['a1'].append(a1)
    res = pd.DataFrame(o)
    res.drop_duplicates(subset=['rsid'], keep=False, inplace=True)
    res['is_ambiguous'] = is_ambiguous(res.a0.tolist(), res.a1.tolist())
    return res

def is_ambiguous(a0, a1):
    o = []
    for i, j in zip(a0, a1):
        if i == PAIRED_BASE[j]:
            o.append(True)
        else:
            o.append(False)
    return o


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='ukb_snp_qc.py', description='''
        Return the list of input SNPs that are non-ambiguous and bi-alleleic.
    ''')
    parser.add_argument('--snplist', help='''
        A list of SNP rsIDs.
    ''')
    parser.add_argument('--ukb-bgi', help='''
        UKB genotype BGEN BGI (from bgenix).
        Need {chr_num} as wildcard.
    ''')
    parser.add_argument('--output', help='''
        Path to the output list.
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
    
    logging.info('Loading SNP list.')
    snplist = read_list(args.snplist)
    snpset = set(snplist)
    
    logging.info('Going through all chromosomes.')
    res = []
    for i in range(1, 23):
        logging.info(f'Working on chromosome{i}.')
        bgi = args.ukb_bgi.format(chr_num=i)
        df_var = extract_var_df(bgi, snpset)
        n0 = df_var.shape[0]
        df_var = df_var[ ~ df_var.is_ambiguous ].reset_index(drop=True)
        logging.info('Working on chromosome {}: before {} SNPs; after {} SNPs'.format(i, n0, df_var.shape[0]))
        res += df_var.rsid.tolist()
    
    logging.info('Writing results.')    
    write_list(args.output, res)
    
