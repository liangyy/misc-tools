This module implement two meta-analysis schemes in R

1. Sample size based
2. Inverse variance based

See [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922887/) Table 1 for more details

Example run:

```
$ Rscript run_meta_analysis.R \
  -G ada_test.txt \
  -a ukbb_test.txt \
  -b Effect,Beta \
  -s StdErr,SE \
  -t Chr:Position \
  -e SNP,Chr,Pos,EA,NEA,EAF \
  -u `pwd` \
  -j annotated_snpid,annotated_snpid \
  -o test-out.txt
```