Quick guide to get your hands on `tabix`.

# Download

Download link: [http://www.htslib.org/download/](http://www.htslib.org/download/) and look for `htslib`.
Look at `INSTALL` note to see how to compile.
At the step of `make install`, you may want to set customized directory for installation by `make install prefix=some-path-here`.

# Quick example

Now that we will rely on some example data in `../annotate_snp_by_position` to see how it works

**First of all**, `tabix` needs files to be compressed with `bgz` (a block-wise gzip shipped as part of `htslib` as well and it is widely adapted in genomics/genetics) which we can treat it just like GZ file but `tabix` will see the difference. 
So, we have to convert the format first.

```
htslibPATH=/gpfs/data/im-lab/nas40t2/yanyul/softwares/htslib-1.10.2
zcat ../annotate_snp_by_position/test_lookup.gz | $htslibPATH/bgzip > test_lookup.bgz 
```

**Next**, `tabix` needs us to specify the format of the input file.
For the purpose of annotation a file with chromosome and position being the first two columns (TAB-delimited), we can pretend it as VCF.

```
$htslibPATH/tabix test_lookup.bgz -p vcf 
```

Once the command is finished (instant!), we will see an index file with extension `tbi`.

```
ls test_lookup*
```
 
**Finally**, we can make query.

```
$htslibPATH/tabix test_lookup.bgz "chr1:149788730-149788790"
```
 
