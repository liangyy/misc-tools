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
parser.add_argument('--snpid_col', type = int, default = 0, help='''
    column index of snpid (CHR:POS) in --input
''')
parser.add_argument('--chr_col', type = int, default = 0, help='''
    column index of CHR in --input
''')
parser.add_argument('--pos_col', type = int, default = 0, help='''
    column index of POS in --input
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
parser.add_argument('--if_input_has_header', default = 1, type = int, help='''
    set to non-1 if does not have header
''')
args = parser.parse_args()

import gzip, re, os


# first go through lookup table and get save position info for each snp 
# then scan the target file

def my_read(filename):
    filep, filee = os.path.splitext(filename)
    if filee == '.gz':
        return gzip.open(filename, 'rt')
    else:
       return open(filename, 'r')

var_dic = {}

with my_read(args.lookup_table) as f:
    for i in f:
        i = i.strip().split('\t')
        chrm = i[args.lookup_chr_col - 1]
        chrm = re.sub('chr', '', chrm)
        start = i[args.lookup_start_col - 1]
        end = i[args.lookup_end_col - 1]
        newid = i[args.lookup_newid_col - 1]
        if start != end:
            continue
        v = chrm + ':' + start
        if v not in var_dic:
            var_dic[v] = newid

o = gzip.open(args.out_txtgz, 'wt')


with my_read(args.input) as f:
    if args.if_input_has_header == 1:
        o.write(next(f).strip() + '\t' + 'new_id' + '\n')
    for i in f:
        i = i.strip().split('\t')
        if args.snpid_col != 0:
            snpid = i[args.snpid_col - 1]
        elif args.chr_col != 0 and args.pos_col != 0:
            snpid = re.sub('chr', '', i[args.chr_col - 1]) + ':' + i[args.pos_col - 1]
        else:
            print("Both --snpid_col and --chr_col/--pos_col are empty! At least one of them should be specified")
            os.exit()
        if snpid in var_dic:
            i.append(var_dic[snpid])
            o.write('\t'.join(i) + '\n')
o.close()



