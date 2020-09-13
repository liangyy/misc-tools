Here we perform LD clumping using plink 1.9.

The input includes:
* genotype (in plink format) 
* GWAS (include p-value and SNP ID that is used in the genotype)

Test run:

```
snakemake -s run.snmk --configfile config.test.yaml -p --config gwas_tag=ADIPOGen_Adiponectin
```
