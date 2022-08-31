import os
import pandas as pd
import argparse
from download_from_ena import ena_ftp

from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool


parser = argparse.ArgumentParser(description='Download data from ENA FTP server')

parser.add_argument('-i','--input', help='File containing sample and reads information', required=True)
parser.add_argument('-o','--output', help='Name of the folder in which all data will be saved', required=True)
parser.add_argument('-n','--threads', type=int, help='Number of threads to use for parallelization, default 1', required=False)
parser.add_argument('-d','--dryrun', help='Dry run, show which runs will be downloaded and their file size', required=False, action='store_true')
#parser.add_argument('-p','--progress', help='Show download progress', required=False, action='store_true')
parser.add_argument('-r','--report_errors', help='Output file report containing failed downloads', required=False)

args = vars(parser.parse_args())

df = pd.read_csv(args['input'])

rows_list = [df.loc[i] for i in df.index]

run_number=len(rows_list)


def download_run(row):
    download_errors = []
    run_data = ena_ftp(project=row['study_accession'], sample=row['accession'], run=row['run_accession'], urls = row['fastq_ftp'], md5=row['fastq_md5'])
    filenames = run_data.files_paths().keys()
    if args['dryrun']:
        print(run_data.run, run_data.files_sizes())
    else:
        for file in filenames:
            downloaded_run = run_data.download_fastq(output_dir=args['output'], filename=file)
            local_md5=run_data.checksum(downloaded_run)
            remote_md5 = run_data.files_md5()[file]

            if local_md5 == remote_md5:
                continue
            else:
                download_errors.append(file)
    return(download_errors)

if args['threads']:
    pool=ThreadPool(args['threads'])
else:
    pool=ThreadPool(1)

download_errors_list = pool.map(download_run, rows_list)

pool.close()
pool.join()

if args['report_errors']:
    download_errors = [j for i in download_errors_list for j in i]
    with open(args['report_errors'],'w') as f:
        for i in download_errors:
            f.write(i)
            f.write('\n')
