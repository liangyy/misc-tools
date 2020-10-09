In the context of UK Biobank genotype data (usually imputed genotype),
for the purpose of some analysis, we may want to limit to SNPs that are:

* Bi-alleleic.
* Non-ambiguous.

Here, I provide a short script to do so.

**Input**:

* A list of SNP rsID.
* UKB BGEN BGI file (it contains the meta data for the variants).

**Output**:

* The list of input SNPs that passed QC.

 