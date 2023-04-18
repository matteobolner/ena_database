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
