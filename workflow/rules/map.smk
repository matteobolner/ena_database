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
    benchmark:
        "benchmarks/bwa_mem2/{sample}/{unit}.txt"
    threads: 4
    wrapper:
        "v1.25.0/bio/bwa-mem2/mem"

rule index_bam:
    input:
        "/lustre/home/bolner/ENA/data/bam/{sample}/{unit}.bam",
    output:
        "/lustre/home/bolner/ENA/data/bam/{sample}/{unit}.bam.bai",
    shell:
        "samtools index {input}"

rule mosdepth:
    input:
        bam="/lustre/home/bolner/ENA/data/bam/{sample}/{unit}.bam"
    output:
        "/lustre/home/bolner/ENA/stats/mosdepth/{sample}/{unit}.mosdepth.summary.txt"
    threads:
        4
    params:
        prefix="/lustre/home/bolner/ENA/stats/mosdepth/{sample}/{unit}"
    shell:
        "mosdepth --fast-mode --no-per-base --threads {threads}  {params.prefix} {input.bam}"

rule samtools_flagstats:
    input:
        bam="/lustre/home/bolner/ENA/data/bam/{sample}/{unit}.bam"
    output:
        "/lustre/home/bolner/ENA/stats/samtools_flagstats/{sample}.tsv"
    threads:
        4
    shell:
        "samtools flagstats -@ {threads} -O tsv {input} > {output}"


rule multiqc_mosdepth:
    input:
        expand("stats/fastp/{u.sample}/{u.unit}_fastp.json",u=units.itertuples())
    output:
        "stats/multiqc/fastp.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/fastp.log"
    wrapper:
        "v1.25.0/bio/multiqc"
