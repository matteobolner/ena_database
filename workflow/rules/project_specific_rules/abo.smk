rule link_fq_to_bam_folder:
    input:
        get_sample_preprocessed_reads
    output:
        fq1="/lustre/home/bolner/PROJECTS/ABO/data/paired_fastq/{sample}.end1.fq",
        fq2="/lustre/home/bolner/PROJECTS/ABO/data/paired_fastq/{sample}.end2.fq"
    shell:
        "ln -s {input[0]} {output.fq1} && ln -s {input[1]} {output.fq2}"


rule link_bam_to_data_folder:
    input:
        bam="/lustre/home/bolner/ENA/data/bam/{sample}.bam",
        idx="/lustre/home/bolner/ENA/data/bam/{sample}.bam.bai"
    output:
        bam="/lustrehome/bolner/work/workspace/ABO_pipeline/resources/samples/{sample}.rmdup.bam",
        idx="/lustrehome/bolner/work/workspace/ABO_pipeline/resources/samples/{sample}.rmdup.bam.bai"
    shell:
        "ln -s {input.bam} {output.bam} && ln -s {input.idx} {output.idx}"
