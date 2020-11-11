import pandas as pd
import os.path

def file_exists(fn):
    return os.path.isfile(fn)

def read_table(fn, indiv_col):
    _, fn_ext = os.path.splitext(fn)
    if fn_ext == '.gz':
        fn_new = re.sub('.gz$', '', fn)
        compress_args = {'compression': 'gzip'}
        _, fn_new = os.path.splitext(fn)
    if fn_ext == '.parquet':
        df = pd.read_parquet(fn)
    elif fn_ext == '.csv':
        df = pd.read_csv(fn, **compress_args)
    elif fn_ext == '.txt' or fn_ext == '.tsv':
        df = pd.read_csv(fn, sep='\s+', **compress_args)
    for i in range(df.shape[1]):
        if df.columns[i] == indiv_col:
            break
    col_list = df.columns.to_list()
    col_list.pop(i)
    col_list = [ indiv_col ] + col_list
    df = df.reindex(columns=col_list)
    df.rename(columns={indiv_col: 'indiv'}, inplace=True)
    df.indiv = df.indiv.astype(str)
    return df

def load_list(fn):
    o = []
    with open(fn, 'r') as f:
        for i in f:
            o.append(i.strip())
    return o

def write_list(mylist, fn):
    with open(fn, 'w') as f:
        for l in mylist:
            f.write(l + '\n')

def intersection(l1, l2):
    a = set(l1)
    a = a.intersection(set(l2))
    return sorted(list(a))
    