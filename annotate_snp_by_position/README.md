Example run:

```
$ python annotate_snp_by_position.py --input test_gwas.gz --snpid_col 1 --lookup_table test_lookup.gz --lookup_chr_col 1 --lookup_start_col 2 --lookup_end_col 3 --lookup_newid_col 4 --out_txtgz test_out.gz --if_input_has_header 0
```

or 

```
$ python annotate_snp_by_position.py --input test_gwas.gz --chr_col 2 --pos_col 3 --lookup_table test_lookup.gz --lookup_chr_col 1 --lookup_start_col 2 --lookup_end_col 3 --lookup_newid_col 4 --out_txtgz test2_out.gz --if_input_has_header 0
```

Example run of `awk_equivalence.sh`

```
$ bash awk_equivalence.sh test_gwas.gz test_lookup.gz
```
