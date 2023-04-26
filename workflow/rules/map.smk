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
        temp("/lustre/home/bolner/ENA/data/bam/{sample}-{unit}.bam"),
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

rule merge_same_sample_bamfiles:
    input:
        get_sample_bams
    output:
        "/lustre/home/bolner/ENA/data/bam/{sample}.bam",
    params:
        bams_per_sample=get_number_of_bams_per_sample
    shell:
        "if [ {params.bams_per_sample} -gt 1 ] ; then samtools merge {input} -o {output}; else cp {input} {output}; fi"

rule index_bam:
    input:
        "/lustre/home/bolner/ENA/data/bam/{sample}.bam",
    output:
        "/lustre/home/bolner/ENA/data/bam/{sample}.bam.bai",
    shell:
        "samtools index {input}"
