import pandas as pd

def union(l1, l2):
    s1 = set(l1)
    return list(s1.union(set(l2)))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='extract.py', description='''
        Extract a list of individuals from a table based on one column.
    ''')
    parser.add_argument('--input', help='''
        Input table TSV with header.
    ''')
    parser.add_argument('--colname', help='''
        column names of the column on which the extraction is based.
    ''')
    parser.add_argument('--by', help='''
        The value to match with.
    ''')
    parser.add_argument('--output_cols', nargs='+', help='''
        Columns to output.
    ''')
    parser.add_argument('--mode', default='match', help='''
        The mode (match/maf) determines how to select the rows.
        Match means to match exactly.
        Maf requires the column to be numeric and will select rows
        with min(value, 1 - value) greater than the by value.
    ''')
    parser.add_argument('--output', help='''
        Output list.
    ''')
    args = parser.parse_args()
 
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    
    logging.info('Loading table.')
    df = pd.read_csv(args.input, sep='\t', header=0)
    pop = args.colname
    if pop not in df.columns or indiv not in df.columns:
        raise ValueError('Colnames are not in input table.')
    colnames = [pop]
    colnames = union(colnames, args.output_cols)
    outnames = args.output_cols
    
    df = df[ colnames ].copy()
    
    logging.info('Extracting.')
    if args.mode == 'match':
        df = df[ df[pop] == args.by ].reset_index(drop=True)
    elif args.mode == 'gt':
        maf = df[pop]
        maf[maf > 0.5] = 1 - df[pop][maf > 0.5]
        df = df[ maf > float(args.gt) ].reset_index(drop=True)
    logging.info('Extracted {} instances.'.format(df.shape[0]))
    
    logging.info('Output.')
    with open(args.output, 'w') as f:
        for i in range(df.shape[0]):
            o = []
            for cc in outnames:
                o.append(df[cc][i])
            f.write('\t'.join(o) + '\n')
    
    logging.info('Done.')
    
    