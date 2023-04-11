import pandas as pd

bs=pd.read_table("data/suina/clean_bs_db.tsv")
df=pd.read_table("data/suina/ena_db.tsv")

units=pd.read_table("data/suina/samples_units.tsv")
units

df=df[~df['scientific_name'].str.contains("scrofa")]
pafr=df[df['scientific_name'].str.contains("Phacochoerus")]
pafr['depth']=pafr['base_count']/2425573090
pafr[(pafr['depth']>5)&(pafr['depth']<25)]['study_accession'].unique()
pafr=pafr[(pafr['study_accession']=='PRJNA837362')&(pafr['depth']>10)]

units=units[units['run'].isin(pafr['run_accession'])]
units.to_csv("config/units.tsv", index=False, sep='\t')

len(pafr)

pafr.to_csv("data/suina/pafr_test.csv")
pafr

df[df['']]
df['scientific_name'].value_counts()
df['scientific_name'].value_counts()


df[df['scientific_name']=='Phacochoerus africanus']
