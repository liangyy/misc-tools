import argparse
parser = argparse.ArgumentParser(prog='liftover_snp.py', description='''
    Lift over from one genome build to another
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
parser.add_argument('--input_delim', type = str, default = 'tab', help='''
    separator in input, space or tab
''')
parser.add_argument('--if_with_chr', default = 0, type = int, help='''
    if you want to add chr before chromosome, e.g. chr1 instead of 1, then set it to 1
''')
args = parser.parse_args()

import pandas as pd
from lib import liftover
import re, os

sep = '\t'
if args.input_delim == 'space':
    sep = ' '
elif args.input_delim == 'tab':
    sep = '\t'
else:
    sep = args.input_delim

filename, file_extension = os.path.splitext(args.input)
if file_extension != 'gz':
    df = pd.read_table(args.input, sep = sep)
else:
    df = pd.read_table(args.input, sep = sep, compression = 'gzip')

chr = df.iloc[:, args.chr_col - 1]
pos = df.iloc[:, args.pos_col - 1] 
# print(chr)
chr = chr.apply(lambda x: re.sub('chr', '', str(x)))
out = liftover(chr, pos, args.liftover_chain)
df = pd.concat([df.reset_index(drop=True), out], axis=1)
if args.if_with_chr != 1:
    df['liftover_chr'] = df['liftover_chr'].apply(lambda x: re.sub('chr', '', x))
df.to_csv(args.out_txtgz, index = False, sep='\t', compression = 'gzip')

