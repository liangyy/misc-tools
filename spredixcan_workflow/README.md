This tends to be an easy S-PrediXcan work pipe.
Everything is derived from [this tutorial in MetaXcan wiki](https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS).

The workflow allows both harmonization and dirty/quick harmonization.


TODO: need to fix the following error:
```
Traceback (most recent call last):
  File "/gpfs/data/im-lab/nas40t2/yanyul/GitHub/MetaXcan/software/M03_betas.py", line 179, in <module>
    run(args)
  File "/gpfs/data/im-lab/nas40t2/yanyul/GitHub/MetaXcan/software/M03_betas.py", line 123, in run
    b = build_betas(args, model, gwas_format, name, args.snp_map_file)
  File "/gpfs/data/im-lab/nas40t2/yanyul/GitHub/MetaXcan/software/M03_betas.py", line 46, in build_betas
    snp_map_ = snp_map.rename(columns={"a0":PF.K_NON_EFFECT_ALLELE, "a1":PF.K_EFFECT_ALLELE})[[PF.K_RSID, PF.K_EFFECT_ALLELE, PF.K_NON_EFFECT_ALLELE, "panel_variant_id", "panel_variant_a0", "panel_variant_a1", "swap"]].drop_duplicates()
  File "/home/t.cri.yliang/miniconda2/envs/metaxcan/lib/python3.7/site-packages/pandas/core/frame.py", line 2806, in __getitem__
    indexer = self.loc._get_listlike_indexer(key, axis=1, raise_missing=True)[1]
  File "/home/t.cri.yliang/miniconda2/envs/metaxcan/lib/python3.7/site-packages/pandas/core/indexing.py", line 1553, in _get_listlike_indexer
    keyarr, indexer, o._get_axis_number(axis), raise_missing=raise_missing
  File "/home/t.cri.yliang/miniconda2/envs/metaxcan/lib/python3.7/site-packages/pandas/core/indexing.py", line 1640, in _validate_read_indexer
    raise KeyError(f"None of [{key}] are in the [{axis_name}]")
KeyError: "None of [Index(['rsid', 'effect_allele', 'non_effect_allele', 'panel_variant_id',\n       'panel_variant_a0', 'panel_variant_a1', 'swap'],\n      dtype='object')] are in the [columns]"
Error in job harmonize while creating output file output/my_gwas_test.harmonization_simple.txt.gz.
RuleException:
CalledProcessError in line 44 of /gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/spredixcan_workflow/spredixcan_pipe.snmk:
Command 'python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/MetaXcan/software/M03_betas.py                 --snp_map_file /gpfs/data/im-lab/nas40t2/abarbeira/projects/gtex_v8/data_formatting/dbsnp/results/snp150_hg19_parsed.txt.gz                 --gwas_file /scratch/t.cri.yliang/predixcan_simulation/test_gwas.txt.gz                 --snp_column SNP                 --non_effect_allele_column REF                 --effect_allele_column ALT                 --beta_column all_inv_var_meta_beta                 --pvalue_column all_inv_var_meta_p                 --keep_non_rsid                 --throw                 --output output/my_gwas_test.harmonization_simple.txt.gz' returned non-zero exit status 1.
  File "/gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/spredixcan_workflow/spredixcan_pipe.snmk", line 44, in __rule_harmonize
  File "/home/t.cri.yliang/miniconda2/envs/mixqtl/lib/python3.6/concurrent/futures/thread.py", line 56, in run
Will exit after finishing currently running jobs.
Exiting because a job execution failed. Look above for error message
``` 
