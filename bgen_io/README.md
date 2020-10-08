Here I build a Python wrapper for UKB BGEN IO. 
The backend is [`bgen-reader`](https://bgen-reader.readthedocs.io/en/latest/index.html).
**Importantly**, we use `bgen-reader v4.0.4` Numpy API with installation as suggested in this [issue](https://github.com/limix/bgen-reader-py/issues/30).
As for 10/8/2020, this wrapper works for un-phased genotype (imputed genotype v3). 
For the phased one, I have one wrapper [here](https://github.com/liangyy/haplotype-po/blob/master/scripts/prs/ukb_hap_reader.py).

# TODO

Migrate the phased BGEN IO to here as well.
