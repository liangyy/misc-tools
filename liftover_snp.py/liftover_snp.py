import argparse
parser = argparse.ArgumentParser(prog='map2gtexv8.py', description='''
    Lift over from hg37 to hg38 
''')
parser.add_argument('--input', help='''
    input file name
''')
parser.add_argument('--chr_col', type = int, help='''
    column index of chromosome (start from 1)
''')
parser.add_argument('--pos_col', type = int, help='''
    column index of position (start from 1)
''')
parser.add_argument('--liftover_chain', help='''
    liftover chain file: e.g. hg17ToHg18.over.chain.gz
''')
parser.add_argument('--out_txtgz', help='''
    output txt.gz file name
''')
args = parser.parse_args()

import pandas as pd
from lib import liftover


df = pd.read_table(args.input)

chr = df.iloc[:, args.chr_col - 1]
pos = df.iloc[:, args.pos_col - 1] 
out = liftover(chr, pos, args.liftover_chain)
df = pd.concat([df.reset_index(drop=True), out], axis=1)
df.to_csv(args.out_txtgz, index = False, sep='\t', compression = 'gzip')

