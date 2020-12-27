import argparse
parser = argparse.ArgumentParser(prog='extract_one_gene.py', description='''
    Extract a list of gene from h5.
''')
parser.add_argument('--input', help='''
    input h5
''')
parser.add_argument('--gene_list', help='''
    path to gene list
''')
parser.add_argument('--output', default = 0, help='''
    output csv
''')
args = parser.parse_args()

import h5py
import numpy as np
import pandas as pd
from pyutil import load_list

gene_list = load_list(args.gene_list)
with h5py.File(args.input, 'r') as f:
    genes = f['genes'][:].astype(str)
    gene_idxs = []
    gene_list_new = []
    for gene_id in gene_list:
        gene_idx = np.where(np.isin(genes, gene_id))[0]
        if len(gene_idx) == 0:
            continue
        else:
            gene_idxs.append(gene_idx)
            gene_list_new.append(gene_id)
    gene_idxs = np.concatenate(gene_idxs, axis=0)
    gene_list = gene_list_new

df_tmp = pd.DataFrame({'idx': gene_idxs, 'gene': gene_list})
df_tmp.sort_values(by='idx', inplace=True)
gene_list = list(df_tmp.gene)
gene_idxs = list(df_tmp.idx)

if len(gene_idxs) == 0:
    print('No target gene in the input.')
else:
    with h5py.File(args.input, 'r') as f:
        pred_expr = f['pred_expr'][gene_idxs, :]
        eid = f['samples'][:].astype(str)
    df = pd.DataFrame({'eid': eid})
    tmp = pd.DataFrame(pred_expr.T, columns=gene_list)
    df = pd.concat([df, tmp], axis=1)
    df.to_csv(args.output, compression='gzip', index=False)


