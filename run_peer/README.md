A simple pipeline for running PEER factor analysis on a give matrix.

To run test.

```
snakemake -s peer.snmk --configfile config.yaml -p
```

*Important*: set `if_transpose` to `Yes` if the matrix is gene x individual


Real example is at [link](https://github.com/liangyy/haplotype-po/blob/master/scripts/framingham_detour/calc_peer/)
