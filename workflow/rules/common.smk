import pandas as pd

configfile: "config/config.yaml"


samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index



##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),

##### Helper functions #####

def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_fastq_remote(wildcards):
    """Get fastq remote path of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def get_raw_reads(wildcards):
    """Get raw reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "/lustre/home/bolner/ENA/data/fastq/raw/{sample}/{unit}_{group}.fastq.gz",
            group=[1, 2],
            **wildcards
        )
    # single end sample
    return "/lustre/home/bolner/ENA/data/fastq/raw/{sample}/{unit}.fastq.gz".format(**wildcards)

def get_preprocessed_reads(wildcards):
    """Get preprocessed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "/lustre/home/bolner/ENA/data/fastq/preprocessed/{sample}/{unit}_{group}.fastq.gz",
            group=[1, 2],
            **wildcards
        )
    # single end sample
    return "/lustre/home/bolner/ENA/data/fastq/preprocessed/{sample}/{unit}.fastq.gz".format(**wildcards)

def get_paired_fastp(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "stats/fastp/pe/{u.sample}/{u.unit}_fastp.json",
            **wildcards
        )
    # single end sample
    return "stats/fastp/se/{u.sample}/{u.unit}_fastp.json".format(**wildcards)


#def get_raw_reads(wildcards):
#    """Get fastq path of given sample-unit."""
#    raw_fastq_path=f"{config['raw_fastq_dir']}/{wildcards.sample}/{wildcards.unit}"
#    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
#    if len(fastqs) == 2:
#        return {"r1": f'{raw_fastq_path}_1.fastq.gz', "r2": f'{raw_fastq_path}_2.fastq.gz'}
#    return {"r1": f'{raw_fastq_path}.fastq.gz'}
#
