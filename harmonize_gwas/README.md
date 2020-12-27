Here we provide a pipe to clean up the GWAS datasets.

The goal is to format these GWAS datasets such that they have columns: 

* variant_id
* effect_allele
* non_effect_allele
* effect size columns:
    1. Set A:
        * effect_size
        * effect_se
    2. Set B:
        * zscore
        * frequency
        * sample_size
* chromosome
* position (optional)


**Features**: 

* Need to handle 'OR'. (x done inside `gwas_parsing`)
* Maybe add AF columns from reference panel. (x done by using reference SNP metadata)
* Handle 'sample_size_half'. (x done inside `gwas_parsing` by assign `n_cases = n_controls = sample_size_half`)

**Procedure**:

1. Use `summary-gwas-imputation/src/gwas_parsing.py` [link](https://github.com/hakyimlab/summary-gwas-imputation/blob/master/src/gwas_parsing.py)
2. Post-process to handle 'OR' columns and 'sample_size_half'


