Here we focus on estimating polygenicity of a complex trait using GWAS summary statistics.

Approaches:

* GCTB - SbayesS [ref](https://www.nature.com/articles/s41467-021-21446-3)
* LD4M [ref](https://www.sciencedirect.com/science/article/pii/S0002929719302666)
* FMR [ref](https://www.biorxiv.org/content/10.1101/2020.09.19.304097v3) Maybe not now ..
 
# SLD4M running mode

* `base`: no stratification on SNPs.
* `all`: include all annotations.
* `all_w_agg`: include all annotations plus an aggregated category by aggregating MAF bins for common variants. 

# Output of SLD4M CSV table

* `*_est`: point estimate, `*_err`: standard error
* `Ma`: polygenicity Me
* `h2`: heritability estimate (maybe scaled?)
* `Maenrich`: estimated Me enrichment (relative to base annotation, i.e. all SNP average)
* `h2enrich`: estimated heritability enrichment (relative to base annotation, i.e. all SNP average)


