gcta_dir=/Users/yanyul/Documents/softwares/gcta_1.92.4beta_mac
pheno_table=test_pheno
grm=test_grm
out_prefix=out

# $gcta_dir/bin/gcta64 --bfile $gcta_dir/test --make-grm-gz --out $grm

# python prep_test.py $gcta_dir/test.phen $pheno_table

for i in `seq 1 4`
do
  echo $gcta_dir/bin/gcta64 --grm-gz $grm --reml \
    --pheno $pheno_table.mphen \
    --mpheno $i \
    --out ${out_prefix}_$i \
    --reml-alg 1
done

for i in `seq 1 4`
do
  echo $gcta_dir/bin/gcta64 --grm-gz $grm --HEreg \
    --pheno $pheno_table.mphen \
    --mpheno $i \
    --out ${out_prefix}_HE_$i
done

export PYTHONPATH=/Users/yanyul/Documents/repo/github/misc-tools/pyutil
python ../run_pyemma.py \
  --grm $grm \
  --grm_cache $grm.cache \
  --y_table $pheno_table.parquet iid \
  --output ${out_prefix}.reml.tsv.gz \
  --reml
  
python ../run_pyemma.py \
  --grm $grm \
  --grm_cache $grm.cache \
  --y_table $pheno_table.parquet iid \
  --output ${out_prefix}.mle.tsv.gz \
    
  