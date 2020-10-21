import sys
import pandas as pd
import numpy as np



infile = sys.argv[1]
outprefix = sys.argv[2]
df = pd.read_csv(infile, sep='\s+', header=None)
df.columns = ['fam', 'iid', 'pheno_1']
df['pheno_2'] = df.pheno_1 + np.random.normal(size=(df.shape[0]))
df['pheno_3'] = df.pheno_1 + np.random.normal(size=(df.shape[0]))
df['pheno_4'] = df.pheno_1 + np.random.normal(size=(df.shape[0]))

cols = [ f'pheno_{i}' for i in range(1, 5) ]
for ph in cols:
    tmp = df[ph].values
    tmp = tmp - tmp.mean()
    df[ph] = tmp
df.to_csv(outprefix + '.mphen', index=False, header=None, sep='\t')
df.drop(columns=['fam'], inplace=True)
df.to_parquet(outprefix + '.parquet', index=False)
