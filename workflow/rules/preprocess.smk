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
        "benchmarks/fastp/{sample}/{unit}.txt"
    shell:
        "fastp -i {input[0]} -I {input[1]} -o {output.fq1} -O {output.fq2} -h {output.html} -j {output.json} -w {threads}"

#rule fastp_se:
#    input:
#        get_raw_reads
#    output:
#        html="stats/fastp/{sample}/{unit}_fastp.html",
#        json="stats/fastp/{sample}/{unit}_fastp.json",
#        fq="/lustre/home/bolner/ENA/data/fastq/preprocessed/{sample}/{unit}.fastq.gz",
#    threads:
#        4
#    benchmark:
#        "benchmarks/fastp/{sample}/{unit}.txt"
#    shell:
#        "fastp -i {input} -o {output.fq} -h {output.html} -j {output.json} -w {threads}"
