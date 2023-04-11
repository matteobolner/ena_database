import pandas as pd
#cd ../
#cd ena_database/

import numpy as np


bs=pd.read_table("data/suina/clean_bs_db.tsv")
df=pd.read_table("data/suina/ena_db.tsv")

df['scientific_name'].value_counts()

df['depth']=df['base_count']/2425573090

porcula=df[df['scientific_name']=='Porcula salvania']

porcula['study_accession']


porcula=porcula[(porcula['depth']>10)&(porcula['depth']<=20)].sort_values(by='depth')
porcula


porcula.to_csv("data/suina/porcula_wur.csv", index=False)



phac=df[df['scientific_name'].str.contains("Phacochoerus")]

phac[phac['depth']>1]

phac['depth']
phac.sort_values(by='depth')
phac[~phac['depth'].isna()]

phac[phac['accession']=='SAMN24148833'].transpose().dropna()
