import argparse
parser = argparse.ArgumentParser(prog='map2gtexv8.py', description='''
    Map variant to SNP ID system given in the lookup table
    (by position and filter out ones with confilct reference/alternative alleles)
    It will report the number of :
      1. exact map;
      2. filpped strand;
      3. filpped ref/alt;
      4. filpped strand and ref/alt
    CAUTION: it only works for SNV but NOT indel!
''')
parser.add_argument('--input', help='''
    input variant list (it should have header!)
''')
parser.add_argument('--chr_col', type = int, help='''
    column index of chromosome (start from 1) in --input
''')
parser.add_argument('--pos_col', type = int, help='''
    column index of position (start from 1) in --input
''')
parser.add_argument('--ref_col', type = int, help='''
    column index of reference allele (start from 1) in --input
''')
parser.add_argument('--alt_col', type = int, help='''
    column index of alternative allele (start from 1) in --input
''')
parser.add_argument('--effect_col', type = str, help='''
    column index of effect (start from 1) in --input which may need to flip the sign.
    separated by ','
''')
parser.add_argument('--lookup_table', help='''
    SNP system lookup table, 
    It will be interpreted as by chromosome 
    if it used something like path/chr@.gz  
''')
parser.add_argument('--lookup_chr_col', type = int, help='''
    column index of chromosome (start from 1) in --lookup_table
''')
parser.add_argument('--lookup_pos_col', type = int, help='''
    column index of position (start from 1) in --lookup_table
''')
parser.add_argument('--lookup_ref_col', type = int, help='''
    column index of reference allele (start from 1) in --lookup_table
''')
parser.add_argument('--lookup_alt_col', type = int, help='''
    column index of alternative allele (start from 1) in --lookup_table
''')
parser.add_argument('--lookup_snpid_col', type = int, help='''
    column index of SNP ID (start from 1) in --lookup_table
''')
parser.add_argument('--out_txtgz', help='''
    output txt.gz file name
''')
parser.add_argument('--include_ambiguious', action = 'store_true', help='''
    output txt.gz file name
''')
args = parser.parse_args()

import gzip
from lib import ref_alt_to_code
from lib import my_read
from lib import read_in_file_name_allowing_by_chromosome


import gzip


# first go through lookup table and get save position info for each snp with ref/alt info (see ref_alt_to_code in lib.py for detail)

var_dic = {}

for lookup_file in read_in_file_name_allowing_by_chromosome(args.lookup_table):
    with my_read(lookup_file) as f:
        next(f)
        for i in f:
            i = i.strip().split('\t')
            chr = i[args.lookup_chr_col - 1]
            start = i[args.lookup_pos_col - 1]
            ref = i[args.lookup_ref_col - 1]
            alt = i[args.lookup_alt_col - 1]
            if len(ref) > 1 or len(alt) > 1:
                continue
            ref_alt_code = ref_alt_to_code(ref, alt)
            v = '*'.join([chr, start, str(ref_alt_code)])
            # print(ref_alt_code, ref, alt)
            v_ambiguious = '*'.join([chr, start, str(-1 * ref_alt_code)])
            if v_ambiguious in var_dic or v in var_dic:
                print('Wrong lookup: same snp has already been included: ', '\t'.join(i))
                continue
            var_dic[v] = ( i[args.lookup_snpid_col - 1], ref, alt )
    print('Finished scanning lookup table: ', lookup_file, flush = True)

meta_dic = {
    'other invalid (long ref or alt)': 0,
    'not shown': 0,
    'exact match': 0,
    'flipped strand': 0,
    'flipped ref/alt': 0,
    'flipped strand and ref/alt': 0,
    'Ambiguious SNP (T/A or G/C SNP)': 0
}
# print(meta_dic)

flip_cols = [ int(i) for i in args.effect_col.split(',') ]

o = gzip.open(args.out_txtgz, 'wt')
with my_read(args.input) as f:
    o.write(next(f).strip() + '\t' + 'annotated_snpid' + '\n')
    for i in f:
        i = i.strip().split('\t')
        chr = i[args.chr_col - 1]
        start = i[args.pos_col - 1]
        ref = i[args.ref_col - 1]
        alt = i[args.alt_col - 1]
        if len(ref) > 1 or len(alt) > 1:
            meta_dic['other invalid (long ref or alt)'] += 1
            continue
        ref_alt_code = ref_alt_to_code(ref, alt)
        if ref_alt_code is 1 or ref_alt_code is 4:
            meta_dic['Ambiguious SNP (T/A or G/C SNP)'] += 1
            if args.include_ambiguious is False:
                continue
        v = '*'.join([chr, start, str(ref_alt_code)])
        v_ambiguious = '*'.join([chr, start, str(-1 * ref_alt_code)])
        if v in var_dic:
            info = var_dic[v]
            newid = info[0]
            if info[1] == ref:
                meta_dic['exact match'] += 1
            else:
                meta_dic['flipped strand'] += 1
                i[args.ref_col - 1] = info[1]
                i[args.alt_col - 1] = info[2]
            i.append(newid)
            o.write('\t'.join(i) + '\n')
        elif v_ambiguious in var_dic:
            info = var_dic[v_ambiguious]
            newid = info[0]
            if info[2] == ref:
                meta_dic['flipped ref/alt'] += 1
            else:
                meta_dic['flipped strand and ref/alt'] += 1
            i.append(newid)
            i[args.ref_col - 1] = info[1]
            i[args.alt_col - 1] = info[2]
            for ff in flip_cols:
                i[ff - 1] = str(-1 * float(i[ff - 1]))
            o.write('\t'.join(i) + '\n')
        else:
            meta_dic['not shown'] += 1
o.close()

for i in list(meta_dic.keys()):
    print('{}\t{}'.format(i, meta_dic[i]))
