import argparse
parser = argparse.ArgumentParser(prog='annotate_snp_by_position.py', description='''
    Map variant with name like CHR:POS to some other ID system by using a lookup table.
    Example data: 
      1. SNPs in input file take the format: chr:pos e.g. 1:100 (with header, TXT.GZ)
      2. Lookup table takes the format: chr start end name, e.g. chr1 10 10 rs1 (TAB separated with/without header in GZ)
    NOTE that the script only consider SNV so that it will ignore all rows in lookup table with start != end. If multiple rows in lookup table have the same position, the first one will be used.
''')
parser.add_argument('--input', help='''
    input variant list 
''')
parser.add_argument('--snpid_col', type = int, help='''
    column index of snpid in --input
''')
parser.add_argument('--lookup_table', help='''
    GTEx V8 variant lookup table 
''')
parser.add_argument('--lookup_chr_col', type = int, help='''
    column index of chromosome (start from 1) in --lookup_table
''')
parser.add_argument('--lookup_start_col', type = int, help='''
    column index of start (start from 1) in --lookup_table
''')
parser.add_argument('--lookup_end_col', type = int, help='''
    column index of reference allele (start from 1) in --lookup_table
''')
parser.add_argument('--lookup_newid_col', type = int, help='''
    column index of new ID (start from 1) in --lookup_table
''')
parser.add_argument('--out_txtgz', help='''
    output txt.gz file name
''')
args = parser.parse_args()

import gzip, re


# first go through lookup table and get save position info for each snp 
# then scan the target file

var_dic = {}

with gzip.open(args.lookup_table, 'rt') as f:
    for i in f:
        i = i.strip().split('\t')
        snpid = i[args.lookup_chr_col - 1]
        chr = re.sub('chr', '', chr)
        start = i[args.lookup_start_col - 1]
        end = i[args.lookup_end_col - 1]
        if start != end:
            continue
        v = chr + ':' + start
        if v not in var_dic:
            var_dic[v] = 1

o = gzip.open(args.out_txtgz, 'wt')


with gzip.open(args.input, 'rt') as f:
    o.write(next(f).strip() + '\t' + 'new_id' + '\n')
    for i in f:
        i = i.strip().split('\t')
        snpid = i[args.snpid_col - 1]
        if snpid in var_dic:
            i.append(var_dic[snpid])
            o.write('\t'.join(i) + '\n')
o.close()



