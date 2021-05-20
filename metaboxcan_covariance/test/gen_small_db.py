import sqlite3
import pandas as pd
db = '/gpfs/data/im-lab/nas40t2/festus/metabolomics/guardian/final_qc/guardian_imp-hapmap-snps_lasso_log-inversenorm_0.01.db'
out = 'small_db.db'
genes = ['V10', 'V101']
conn = sqlite3.connect(db)
df = pd.read_sql_query("SELECT * FROM weights", conn)
df = df[ df.gene.isin(genes) ].reset_index(drop=True)
df2 = pd.read_sql_query("SELECT * FROM extra", conn)
df2 = df2[ df2.gene.isin(genes) ].reset_index(drop=True)
df.to_csv('small_db.tsv', sep='\t')

conn2 = sqlite3.connect(out)
df.to_sql(name='weights', con=conn2)
df2.to_sql(name='extra', con=conn2)
conn.close()
conn2.close()