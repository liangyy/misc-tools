Liftover SNP (implemented by `pyliftover`)

Example run

```
$ python liftover_snp.py --input test-gwas.txt --chr_col 2 --pos_col 3 --liftover_chain ~/labshare/data/hg18ToHg19.over.chain.gz --input_delim space --out_txtgz test-out2.gz --if_with_chr 1
$ python liftover_snp.py --input test-gwas.txt --chr_col 2 --pos_col 3 --liftover_chain ~/labshare/data/hg18ToHg19.over.chain.gz --input_delim space --out_txtgz test-out.gz
```
