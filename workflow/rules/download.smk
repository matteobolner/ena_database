#from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
#FTP = FTPRemoteProvider()

rule download_fq_se:
    params:
        fqs=get_fastq_remote
    output:
        "/lustre/home/bolner/ENA/data/fastq/raw/{sample}_{unit}.fastq.gz"
    benchmark:
        "stats/benchmarks/download/{sample}-{unit}.txt"
    log:
        "logs/download/wget/{sample}-{unit}.log",
    threads:
        1
    shell:
        "wget {params.fqs[r1]} -o {log} -O {output}"

rule download_fq_pe:
    params:
        fqs=get_fastq_remote
    output:
        fq1="/lustre/home/bolner/ENA/data/fastq/raw/{sample}-{unit}.1.fastq.gz",
        fq2="/lustre/home/bolner/ENA/data/fastq/raw/{sample}-{unit}.2.fastq.gz",
    benchmark:
        "stats/benchmarks/download/{sample}-{unit}.txt"
    log:
        fq1="logs/download/wget/{sample}-{unit}.1.log",
        fq2="logs/download/wget/{sample}-{unit}.2.log"
    threads:
        1
    shell:
        "wget {params.fqs[r1]} -o {log.fq1} -O {output.fq1} && wget {params.fqs[r2]} -o {log.fq2} -O {output.fq2}"


rule fastp:
    input:
        fq1="/lustre/home/bolner/ENA/data/fastq/raw/{sample}-{unit}.1.fastq.gz",
        fq2="/lustre/home/bolner/ENA/data/fastq/raw/{sample}-{unit}.2.fastq.gz",
    output:
        html="stats/fastp/{sample}-{unit}_fastp.html",
        json="stats/fastp/{sample}-{unit}_fastp.json",
        fq1="/lustre/home/bolner/ENA/data/fastq/preprocessed/{sample}-{unit}.1.fastq.gz",
        fq2="/lustre/home/bolner/ENA/data/fastq/preprocessed/{sample}-{unit}.2.fastq.gz",
    threads:
        4
    benchmark:
        "stats/benchmarks/fastp/{sample}-{unit}.txt"
    shell:
        "fastp -i {input.fq_1} -I {input.fq_2} -o {output.fq_1} -O {output.fq_2} -h {output.html} -j {output.json} -w {threads}"



rule test_list:
    input:
        reads=get_raw_reads
    output:
        "data/test_list/{sample}-{unit}.txt"
    shell:
        "echo {input} > {output}"
'''


rule split_db_by_sample:
    input:
        db="data/test/test_db.csv"
    output:
        sample_db="data/test_db/{run}.csv"
    run:
        tempdf=pd.read_csv(input.db)
        for name,group in tempdf.groupby(by='run_accession'):
            group.to_csv(f'{output.sample_db}', index=False)


rule test_download:
    input:
        sample="data/test_db/{sample}.csv",
        #sample="data/suina/porcula_wur_4.csv"
    output:
        dir=directory("/lustre/home/bolner/FASTQ/test/{sample}"),
    threads:
        4
    log:
        "logs/download/{sample}.txt"
    benchmark:
        "stats/benchmarks/download/{sample}.txt"
    shell:
        "python3 workflow/scripts/download_from_ena.py -i {input.sample} -o {output.dir} -n {threads} -r {log}"

'''
'''

rule test_wget_speed:
    input:
        "data/test/{run}.csv"
    output:
        "data/test_wget/"



rule test_fastp:
    input:
        fq_1="/lustre/home/bolner/FASTQ/porcula_test/{sample}/{run}_1.fastq.gz",
        fq_2="/lustre/home/bolner/FASTQ/porcula_test/{sample}/{run}_2.fastq.gz"
    output:
        html="stats/fastp/{sample}/{run}_fastp.html",
        json="stats/fastp/{sample}/{run}_fastp.json",
        fq_1="/lustre/home/bolner/FASTP/porcula_test/{sample}/{run}_1.fastq.gz",
        fq_2="/lustre/home/bolner/FASTP/porcula_test/{sample}/{run}_2.fastq.gz",
    threads:
        4
    benchmark:
        "stats/benchmarks/fastp/{sample}_{run}.txt"
    shell:
        "fastp -i {input.fq_1} -I {input.fq_2} -o {output.fq_1} -O {output.fq_2} -h {output.html} -j {output.json} -w {threads}"

rule multiqc_fastp:
    input:
        expand("stats/fastp/{sample}/{run}_fastp.json", zip, sample=samples, run=runs)
    output:
        "stats/fastp/fastp_multiqc.html"
    params:
        dir="stats/fastp/"
    shell:
        "multiqc {params.dir} --filename {output}"

rule bwa_index_ref:
    input:
        "{genome}",
    output:
        "{genome}.0123",
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.bwt.2bit.64",
        "{genome}.pac",
    log:
        "logs/bwa-mem2_index/{genome}.log",
    wrapper:
        "v1.25.0/bio/bwa-mem2/index"

rule bwa_mem2_mem:
    input:
        reads=["/lustre/home/bolner/FASTP/porcula_test/{sample}/{run}_1.fastq.gz", "/lustre/home/bolner/FASTP/porcula_test/{sample}/{run}_1.fastq.gz"],
        idx=multiext(config['ref_genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac")
    output:
        "/lustre/home/bolner/ENA/porcula_test/{sample}_{run}.bam",
    log:
        "logs/bwa_mem2/{sample}_{run}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:unknown\tLB:{sample}'",
        sort="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 4
    wrapper:
        "v1.25.0/bio/bwa-mem2/mem"
'''
