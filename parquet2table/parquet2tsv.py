import argparse
parser = argparse.ArgumentParser(prog='parquet2tsv.py', description='''
    Convert parquet file to TSV
''')
parser.add_argument('--parquet', required=True, help='''
    Input parquet file
''')
parser.add_argument('--output', required=True, help='''
    Path to output file
''')
parser.add_argument('--col_list', nargs='+', default=None, help='''
    A list of columns to extract.
    Add extract column names.
''')
args = parser.parse_args()

import pandas as pd

table = pd.read_parquet(args.parquet)

if args.col_list is None:
    pass
else:
    cols = args.col_list[1:]
    list_cols = []
    with open(args.col_list[0], 'r') as f:
        for i in f:
            i = i.strip()
            list_cols.append(i)
    table = table[ cols + list_cols ]
    
table.to_csv(args.output, sep='\t', index=False)
