Given SNP position and allele information, annotate with SNP ID given in lookup table

It is useful when we want to annotate GWAS summary statistic with the SNP ID system we like to use.

Example run

```
$ python map_by_lookup_table.py --input gwas.txt --chr_col 2 --pos_col 3 --ref_col 4 --alt_col 5 --effect_col 7 --lookup_table lookup.txt --lookup_chr_col 1 --lookup_pos_col 2 --lookup_ref_col 3 --lookup_alt_col 4 --lookup_snpid_col 5 --out_txtgz test.gz
```