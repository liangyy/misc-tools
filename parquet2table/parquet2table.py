import argparse
parser = argparse.ArgumentParser(prog='parquet2table.py', description='''
    Convert parquet file to table
''')
parser.add_argument('--parquet', required=True, help='''
    Input parquet file
''')
parser.add_argument('--output', required=True, help='''
    Path to output file
''')
args = parser.parse_args()

import pyarrow.parquet as pq

table = pq.read_table(args.parquet)
table.to_pandas().to_csv(args.output, compression = 'gzip', sep = '\t', index = False)
