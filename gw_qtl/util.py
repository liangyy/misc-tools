import pandas as pd

from pyutil import read_table, read_yaml


def rearrange_rows(df, target_list):
    tmp = pd.DataFrame({'indiv': target_list})
    return pd.merge(tmp, df, on='indiv', how='left')
    
def exclude_b_from_a(a, b):
    a_ = set(a)
    b_ = set(b)
    a_ = a_.difference(b_)
    return list(a_)

def check_a_in_b(a, b):
    return all(i in b for i in a)

def load_table(file_list):
    fn, indiv_col = file_list[:2]
    return read_table(fn, indiv_col)

def load_covariate(file_list, fyaml):
    df = load_table(file_list)
    df_yaml = read_yaml(fyaml)
    if not check_a_in_b(list(df_yaml.keys()), df.columns.to_list()):
        raise ValueError(f'The {file_string} does not match {fyaml}.')
    df_res = [ pd.DataFrame({'indiv': df.indiv}) ]
    nsample = df.shape[0]
    for col, ctype in df_yaml.items():
        if ctype == 'continuous':
            df_res.append(pd.DataFrame({col: df[col]}))
        if ctype == 'categorical':
            col_values = df[col].values
            categories = np.unique(col_values)
            ncate = categories.shape[0]
            if ncate > 25:
                raise ValueError(f'Too many categories at {col}.')
            mat = np.zeros((nsample, ncate))
            for i in range(ncate):
                mat[ col_values == categories[i], i ] = 1
            df_res.append(
                pd.DataFrame(
                    mat, 
                    columns=[ f'{col}_{val}' for val in categories.tolist() ]
                )
            )
    return pd.concat(df_res, axis=1)

def load_phenotype(file_list, fyaml):
    df = load_table(file_list)
    df = df.reindex(columns=['indiv'] + list(df_yaml.keys()))
    df.drop_duplicates(subset='indiv', inplace=True)
    return df
