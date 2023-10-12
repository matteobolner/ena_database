import pandas as pd
import numpy as np
import seaborn as sns

bs=pd.read_table("data/suina/clean_bs_db.tsv")
df=pd.read_table("data/suina/ena_db.tsv")

df=df[df['library_source']=='GENOMIC']
df=df[df['library_strategy']=='WGS']
df=df[~df['base_count'].isna()]

#-choerus
choerus=df[df['scientific_name'].str.contains("choerus")]
choerus['depth']=choerus['base_count']/2425573090
choerus=choerus[choerus['depth']>5]
choerus=choerus[choerus['depth']<30]

phacochoerus_aethiopicus=choerus[choerus['scientific_name']=='Phacochoerus aethiopicus']

potamochoerus=choerus[choerus['scientific_name'].str.contains("Potamochoerus")]

phacochoerus_africanus=choerus[choerus['scientific_name']=='Phacochoerus africanus']
phacochoerus_africanus=phacochoerus_africanus[phacochoerus_africanus['study_accession']!='PRJNA691462']
phacochoerus_africanus=phacochoerus_africanus.sort_values(by='depth')
phacochoerus_africanus=phacochoerus_africanus.iloc[0:15]

choerus_output=pd.concat([potamochoerus, phacochoerus_aethiopicus, phacochoerus_africanus ])
choerus_output['fq1']=choerus_output['fastq_ftp'].str.split(";").str[0]
choerus_output['fq2']=choerus_output['fastq_ftp'].str.split(";").str[1]
choerus_output=choerus_output[~choerus_output['fastq_ftp'].isna()]
choerus_output.to_csv("data/suina/choerus.csv", index=False)

choerus_output=choerus_output[['accession','run_accession','fq1','fq2']]
choerus_output.columns=['sample','unit','fq1','fq2']
choerus_output.to_csv("config/units.tsv", sep='\t', index=False)
#choerus_output.to_csv("data/suina/choerus.tsv", sep='\t', index=False)
choerus_output[['sample']].to_csv("config/samples.tsv", sep='\t', index=False)


samples=pd.read_table("stats/sample_depths.tsv")

samples=choerus[['scientific_name','sample_accession']].merge(samples, left_on='sample_accession', right_on='sample')
samples=samples.drop(columns=['sample_accession'])

species_name={'Phacochoerus africanus':'African warthog', 'Potamochoerus larvatus':'Bushpig','Potamochoerus porcus':'Red river hog','Phacochoerus aethiopicus':'Desert warthog'}

samples['breed']=samples['scientific_name'].apply(lambda x:species_name[x])
samples['dataset']='ENA'
samples=samples[['sample','breed','dataset','depth','scientific_name']]
samples=samples.rename(columns={'scientific_name':'species'})
samples['origin']="Africa"
samples=samples.sort_values(by='species')
samples.to_csv("data/suina/choerus_sample_info.tsv", index=False, sep='\t')
