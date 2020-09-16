Here we provide a simple Python script to upload a list of files to zenodo.
Internally, it runs zenodo API documented at [here](https://developers.zenodo.org/).

You should have personal token to be ready to go. Here we put the token at `~/token.txt`.

Example:

Print help message.
```
$ python zenodo_uploader.py -h 
usage: zenodo_uploader.py [-h] [--depository DEPOSITORY]
                          [--meta-yaml META_YAML] [--file-list FILE_LIST]
                          [--token TOKEN]

Upload a list of files to zenodo. It will create one zenodo deposite and
upload all files there.

optional arguments:
  -h, --help            show this help message and exit
  --depository DEPOSITORY
                        If want to append to existing depository, set the ID
                        here. Otherwise, new depository will be generated.
  --meta-yaml META_YAML
                        A yaml file containing some meta data to fill along
                        with the files.
  --file-list FILE_LIST
                        The list of files to upload.
  --token TOKEN         Token saved in a file. Note: make sure to delete the
                        token after the uploading is done.
```

Add a list of files to a new depository. The depository ID will be printed out. You can check the uploaded files at URL: https://zenodo.org/deposit/[depository-ID]. 
```
$ python zenodo_uploader.py --meta-yaml meta_info.yaml  --file-list test_files.txt --token ~/zenodo.txt
2020-09-15 10:11:09 PM  Use depository = https://zenodo.org/api/deposit/depositions/[depository-ID]
2020-09-15 10:11:09 PM  Processing mixqtl.Whole_Blood_GTEx_eGene.cis_qtl_pairs.mixQTL.chr20.parquet
```

Add additional files to the depository.
```
$ python zenodo_uploader.py --depository [depository-ID]  --file-list test_files2.txt --token ~/zenodo.txt
2020-09-15 10:12:29 PM  Use depository = https://zenodo.org/api/deposit/depositions/[depository-ID]
2020-09-15 10:12:29 PM  Processing mixqtl.Whole_Blood_GTEx_eGene.cis_qtl_pairs.mixQTL.chr22.parquet
```

Done.

