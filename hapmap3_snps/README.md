This module extract HapMap 3 SNPs for downstream use.
Download data from [https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/]

Dependency: 

* python3: pandas
* plink1.9
* snakemake
* download the repo: [https://github.com/liangyy/misc-tools](https://github.com/liangyy/misc-tools)

Test run:

```
# maf = 0.01
snakemake -s hapmap.snmk --configfile config.hapmap3_eur.yaml -p --config maf=0.01
# maf = 0.05
snakemake -s hapmap.snmk --configfile config.hapmap3_eur.yaml -p --config maf=0.05
```

The resulting files of the test run are here: [maf = 0.01](https://uchicago.box.com/s/4r2ddcyaizmutm0nmfrfg5vr7xz7hsb5) and [maf = 0.05](https://uchicago.box.com/s/7dpcj59y2u7titrevwu76rzudpq2gxfu). 


To generate metadata table for lookup, see `gen_lookup_table.py`.