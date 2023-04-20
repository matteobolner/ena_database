import pandas as pd

bs=pd.read_table("data/suina/clean_bs_db.tsv")
df=pd.read_table("data/suina/ena_db.tsv")

#units=pd.read_table("data/suina/samples_units.tsv")


df=df[~df['scientific_name'].str.contains("scrofa")]

df['scientific_name'].value_counts()

porcula=df[df['scientific_name']=='Porcula salvania']

porcula['depth']=porcula['base_count']/2472047747
porcula['depth']

porcula=porcula[porcula['depth']>10][['sample_accession','run_accession','fastq_ftp']]
porcula['fastq_ftp'].iloc[0]

porcula['fq1']=porcula['fastq_ftp'].str.split(";").str[0]
porcula['fq2']=porcula['fastq_ftp'].str.split(";").str[1]

outdb=df[df['sample_accession'].isin(samples['sample'])]
outdb.to_csv("data/suina/pafr_and_porcula.csv")

outdb=db


######
test=pd.read_csv("/lustrehome/bolner/work/workspace/test/plink_test/test.fam", sep=' ', header=None)

samples=pd.read_csv("config/samples.tsv", sep='\t')
origin={i:j for i,j in zip(samples['sample'], samples['origin'])}
test[0]=test[0].apply(lambda x:origin[x])
test.to_csv("/lustrehome/bolner/work/workspace/test/plink_test/test.fam", header=None, index=False, sep=' ')
######

porcula=porcula[['sample_accession','run_accession','fq1','fq2']]

samples=pd.read_table("config/samples.tsv")

units=pd.read_table("config/units.tsv")

units=units.iloc[0:6]
porcula.columns=units.columns

units=pd.concat([units, porcula])

units.to_csv("config/units.tsv", index=False, sep='\t')
samples=samples.iloc[0:6]
samples=units[['sample']]
samples.to_csv("config/samples.tsv", index=False, sep='\t')

pafr=df[df['scientific_name'].str.contains("Phacochoerus")]
pafr['depth']=pafr['base_count']/2425573090
pafr[(pafr['depth']>5)&(pafr['depth']<25)]['study_accession'].unique()
pafr=pafr[(pafr['study_accession']=='PRJNA837362')&(pafr['depth']>10)]

units=units[units['run'].isin(pafr['run_accession'])]
#units.to_csv("config/units.tsv", index=False, sep='\t')

len(pafr)

pafr.to_csv("data/suina/pafr_test.csv")
pafr

df[df['']]
df['scientific_name'].value_counts()
df['scientific_name'].value_counts()


df[df['scientific_name']=='Phacochoerus africanus']
