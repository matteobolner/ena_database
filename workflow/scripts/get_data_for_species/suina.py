import pandas as pd
import sys
sys.path.append("/lustrehome/bolner/work/workspace/ena_database/scripts/")
from workflow.scripts.build_db_class import ena_db

#db=ena_db(taxon=35497, breeds='resources/breed_names/pig.txt', sexes=['male','female','gilt','sow','gilts','sows','barrow','queen'])

db=ena_db(taxon=35497, breeds='resources/breed_names/pig.txt', sexes=['male','female','gilt','sow','gilts','sows','barrow','queen'])

#db.ena_df.to_csv("data/suina/ena_db.tsv", sep='\t', index=False)
bs_df=db.build_biosamples_df(db.ena_df, thread_number=4)
bs_df.to_csv("data/suina/bs_db.tsv", sep='\t', index=False)

clean_bs=db.clean_biosamples_df(bs_df)
#db.build_biosamples_df(db.ena_df, thread_number=4)
clean_bs.to_csv("data/suina/clean_bs_db.tsv", sep='\t', index=False)


db['fq1']=db['fastq_ftp'].str.split(";").str[0]
db['fq2']=db['fastq_ftp'].str.split(";").str[1]
db=db[~db['fastq_ftp'].isna()]
samples_units=df[['accession','run_accession','fq1','fq2']]
samples_units.columns=['sample','run','fq1','fq2']
samples_units.to_csv("data/suina/samples_units.tsv", sep='\t', index=False)
