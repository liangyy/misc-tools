export PYTHONPATH=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/transethnic_prs/:/gpfs/data/im-lab/nas40t2/yanyul/GitHub/SPrediXcan2PTRS/
conda activate transethnic_prs

genotype=/gpfs/data/im-lab/nas40t2/festus/metabolomics/guardian/final_qc/guardian.imp.qc_filtered.bed
small_db=small_db.db
if [[ ! -f $small_db ]]
then
  python gen_small_db.py 
fi

python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/generate_covariance.py \
  --genotype_bed $genotype \
  --predictdb $small_db \
  --snpid_to_report rsid \
  --output rsid_output.txt.gz

  
python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/generate_covariance.py \
  --genotype_bed $genotype \
  --predictdb $small_db \
  --snpid_to_report rsid varID \
  --output rsid_output