This module extract HapMap 3 SNPs for downstream use.
Download data from [https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3/plink_format/](https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3/plink_format/)

Dependency: 

* python3: see imlabtools conda env at [here](https://github.com/hakyimlab/MetaXcan#example-conda-environment-setup)
* plink1.9
* snakemake
* download the repo: [https://github.com/liangyy/misc-tools](https://github.com/liangyy/misc-tools)

Test run:

```
# build b36
## maf = 0.01
snakemake -s hapmap.snmk --configfile config.hapmap3_eur.yaml -p --config maf=0.01 target_build=b36
## maf = 0.05
snakemake -s hapmap.snmk --configfile config.hapmap3_eur.yaml -p --config maf=0.05 target_build=b36

# build b37
## maf = 0.01
snakemake -s hapmap.snmk --configfile config.hapmap3_eur.yaml -p --config maf=0.01 target_build=b37
## maf = 0.05
snakemake -s hapmap.snmk --configfile config.hapmap3_eur.yaml -p --config maf=0.05 target_build=b37
```

The resulting files of the test run are here: 
* [maf = 0.01, build = b36](https://uchicago.box.com/s/4r2ddcyaizmutm0nmfrfg5vr7xz7hsb5)
* [maf = 0.05, build = b36](https://uchicago.box.com/s/7dpcj59y2u7titrevwu76rzudpq2gxfu)
* [maf = 0.01, build = b37](https://uchicago.box.com/s/4r2ddcyaizmutm0nmfrfg5vr7xz7hsb5)
* [maf = 0.05, build = b37](https://uchicago.box.com/s/7dpcj59y2u7titrevwu76rzudpq2gxfu) 

