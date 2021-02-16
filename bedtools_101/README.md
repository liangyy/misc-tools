Consider a use case where we have a list of LD blocks in `ldblock.bed` and a list of SNPs in `snp.bed` and we want to know which is in which LD block.

```
bedtools intersect -a snp.bed -b ldblock.bed -loj
```

**Explanation**:

* `-a snp.bed` since we treat SNPs as the query (we want to annotate SNPs in some way).
* `-b ldblock.bed` since we could treat the LD block information as the database that we want to query from.
* `-loj` since we also want to keep SNPs that do not overlap with any LD block.

See more information at [here](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)
