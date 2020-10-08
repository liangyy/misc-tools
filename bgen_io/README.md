Here I build a Python wrapper for UKB BGEN IO. 

# Notes on backend

As a python wrapper, I wish to use Python package. 
Initially, I tried setting up the backend as [`bgen-reader`](https://bgen-reader.readthedocs.io/en/latest/index.html).
I first tried Dash based API and it has issues for multithreading mode.
And then, I used `bgen-reader v4.0.4` Numpy API with installation as suggested in this [issue](https://github.com/limix/bgen-reader-py/issues/30).
But it takes too long to initialize.

With such failure, I have to go back to `rbgen`.

# About 

As for 10/8/2020, this wrapper works for un-phased genotype (imputed genotype v3). 
For the phased one, I have one wrapper [here](https://github.com/liangyy/haplotype-po/blob/master/scripts/prs/ukb_hap_reader.py).

# TODO

Migrate the phased BGEN IO to here as well.
