import pandas as pd
from ftplib import FTP
from workflow.scripts.ena import EnaFtp

df=pd.read_csv("data/test/test_db.csv")
df['library_layout']
a=df.transpose()

df.columns

def setup_ftp(root="ftp.sra.ebi.ac.uk", main_dir="vol1/fastq/"):
    ftp = FTP(root)
    ftp.login()
    ftp.cwd(main_dir)
    return ftp



test_ftp=setup_ftp()

for index,row in df.iterrows():

    run_data = EnaFtp(
        project=row["study_accession"],
        sample=row["accession"],
        run=row["run_accession"],
        urls=row["fastq_ftp"],
        md5=row["fastq_md5"],
        test_ftp=test_ftp
    )
    print(row['accession'])
    print(run_data.files_paths())
run_data.urls
