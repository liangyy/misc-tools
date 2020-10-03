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
args = parser.parse_args()

import pandas as pd

table = pd.read_parquet(args.parquet)
table.to_csv(args.output, compression='gzip', sep='\t', index=False)
