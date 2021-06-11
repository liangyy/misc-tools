* Update on 6/11/21: 
    - Add look-up table as the final output.
    - We switch from https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/ to https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3/plink_format/ for the test run. But we still share the resulting files from the old runs: [maf = 0.01](https://uchicago.box.com/s/4r2ddcyaizmutm0nmfrfg5vr7xz7hsb5) and [maf = 0.05](https://uchicago.box.com/s/7dpcj59y2u7titrevwu76rzudpq2gxfu)
    - For record, the old config is

```
name_tag: 'hapmap3_eur'
pop_tag: 'CEU'
individual_list:
  file: '/gpfs/data/im-lab/nas40t2/yanyul/data/hapmap3/relationships_w_pops_121708.txt'
  pop_col: 'population'
  indiv_col: 'IID'
  fam_col: 'FID'
genotypes:
  ped: '/gpfs/data/im-lab/nas40t2/yanyul/data/hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.ped'
  map: '/gpfs/data/im-lab/nas40t2/yanyul/data/hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.map'
plink1.9: 'plink'  # on CRI: module load gcc/6.2.0; module load plink/1.90
plink_mem_in_mb: 10000
misc_tools_path: '/gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/'

gen_lookup_table:
  chain_file: null
  target_build: 'bX'
```