"""
PREPROCESSING
"""
rule multiqc_fastp:
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


"""
MAPPING
"""

rule mosdepth:
    input:
        bam="/lustre/home/bolner/ENA/data/bam/{sample}.bam",
        idx="/lustre/home/bolner/ENA/data/bam/{sample}.bam.bai"
    output:
        "stats/mosdepth/{sample}/{sample}.mosdepth.summary.txt"
    threads:
        4
    params:
        prefix="stats/mosdepth/{sample}/{sample}"
    shell:
        "mosdepth --fast-mode --no-per-base --threads {threads} {params.prefix} {input.bam}"

rule multiqc_mosdepth:
    input:
        expand("stats/mosdepth/{sample}/{sample}.mosdepth.summary.txt", sample=samples.index)
    output:
        "stats/multiqc/mosdepth.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/mosdepth.log"
    wrapper:
        "v1.25.0/bio/multiqc"

rule get_overall_sample_depth:
    input:
        depths=expand("stats/mosdepth/{sample}/{sample}.mosdepth.summary.txt", sample=samples.index),
    output:
        depths="stats/sample_depths.tsv"
    params:
        samples=samples.index
    run:
        outdf=pd.DataFrame(columns=['sample','depth'])
        for i,j in zip(input.depths, params.samples):
            tempdf=pd.read_table(i).set_index('chrom')
            sample_depth=tempdf.loc['total']['mean']
            outdf.loc[len(outdf)]=[j, sample_depth]
        outdf.to_csv(output.depths, sep='\t', index=False)

rule samtools_flagstats:
    input:
        bam="/lustre/home/bolner/ENA/data/bam/{sample}.bam"
    output:
        "stats/samtools_flagstats/{sample}.txt"
    threads:
        4
    shell:
        "samtools flagstats -@ {threads} {input} > {output}"

rule multiqc_flagstats:
    input:
        expand("stats/samtools_flagstats/{sample}.txt", sample=samples.index)
    output:
        "stats/multiqc/samtools_flagstats.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/samtools_flagstats.log"
    wrapper:
        "v1.25.0/bio/multiqc"
