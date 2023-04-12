#from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
#FTP = FTPRemoteProvider()
'''
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
'''

rule split_db_by_unit:
    input:
        #db="data/suina/pafr_test.csv"
        db="data/suina/test.tsv"
    output:
        outfile="data/split_db/{sample}/{unit}.tsv"
    params:
        sample= lambda wc:wc.sample,
        unit= lambda wc:wc.unit
    run:
        tempdf=pd.read_table(input.db)
        tempdf=tempdf[(tempdf['accession']==params.sample)&(tempdf['run_accession']==params.unit)]
        tempdf.to_csv(output.outfile, index=False, sep='\t')

rule download_fastq_se:
    input:
        "data/split_db/{sample}/{unit}.tsv"
    output:
        "/lustre/home/bolner/ENA/data/fastq/raw/{sample}/{unit}.fastq.gz"
    params:
        outdir="/lustre/home/bolner/ENA/data/fastq/raw/"
    threads:
        1
    shell:
        "python3 workflow/scripts/download_from_ena.py -i {input} -o {params.outdir} -n {threads}"

rule download_fastq_pe:
    input:
        "data/split_db/{sample}/{unit}.tsv"
    output:
        "/lustre/home/bolner/ENA/data/fastq/raw/{sample}/{unit}_1.fastq.gz"
        "/lustre/home/bolner/ENA/data/fastq/raw/{sample}/{unit}_2.fastq.gz"
    params:
        outdir="/lustre/home/bolner/ENA/data/fastq/raw/"
    threads:
        1
    shell:
        "python3 workflow/scripts/download_from_ena.py -i {input} -o {params.outdir} -n {threads}"


rule fastp_pe:
    input:
        get_raw_reads
    output:
        html="stats/fastp/{sample}/{unit}_fastp.html",
        json="stats/fastp/{sample}/{unit}_fastp.json",
        fq1="/lustre/home/bolner/ENA/data/fastq/preprocessed/{sample}/{unit}_1.fastq.gz",
        fq2="/lustre/home/bolner/ENA/data/fastq/preprocessed/{sample}/{unit}_2.fastq.gz",
    threads:
        4
    benchmark:
        "stats/benchmarks/fastp/{sample}/{unit}.txt"
    shell:
        "fastp -i {input[0]} -I {input[1]} -o {output.fq1} -O {output.fq2} -h {output.html} -j {output.json} -w {threads}"

rule fastp_se:
    input:
        get_raw_reads
    output:
        html="stats/fastp/{sample}/{unit}_se_fastp.html",
        json="stats/fastp/{sample}/{unit}_se_fastp.json",
        fq="/lustre/home/bolner/ENA/data/fastq/preprocessed/{sample}/{unit}.fastq.gz",
    threads:
        4
    benchmark:
        "stats/benchmarks/fastp/{sample}/{unit}.txt"
    shell:
        "fastp -i {input} -o {output.fq} -h {output.html} -j {output.json} -w {threads}"

rule multiqc_fastp:
    input:
        expand("stats/fastp/{u.sample}/{u.unit}_fastp.json",u=units.itertuples())
    output:
        "stats/multiqc/fastp.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.25.0/bio/multiqc"

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
        reads=get_preprocessed_reads,
        idx=multiext(config['ref_genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac")
    output:
        "/lustre/home/bolner/ENA/data/bam/{sample}/{unit}.bam",
    log:
        "logs/bwa_mem2/{sample}/{unit}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:unknown\tLB:{sample}'",
        sort="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 4
    wrapper:
        "v1.25.0/bio/bwa-mem2/mem"
