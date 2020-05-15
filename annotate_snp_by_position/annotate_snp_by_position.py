import argparse
parser = argparse.ArgumentParser(prog='extract_one_gene.py', description='''
    Extract one gene from h5.
''')
parser.add_argument('--input', help='''
    input h5
''')
parser.add_argument('--gene_id', type = str, help='''
    gene id
''')
parser.add_argument('--output', type = int, default = 0, help='''
    output csv
''')
args = parser.parse_args()

import h5py
import numpy as np
import pandas as pd

with h5py.File(args.input, 'r') as f:
    genes = f['genes'][:].astype(str)
    gene_idx = np.where(np.isin(f, args.gene_id))

if len(gene_idx) == 0:
    print('No {} in the input'.format(args.gene_id))
elif len(gene_idx) == 1:
    gene_idx = gene_idx[0]
    with h5py.File(args.input, 'r') as f:
        pred_expr = f['pred_expr'][gene_idx, :]
        eid = f['samples'][:].astype(str)
    df = pd.DataFrame({'eid': eid, args.gene_id: pred_expr})
    df.to_csv(args.output, compression='gzip')
