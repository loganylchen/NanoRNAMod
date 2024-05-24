rule create_dict:
    input:
        config['reference']['transcriptome_fasta']
    output:
        config['reference']['transcriptome_fasta']+'.dict'
    log:
        "logs/picard/create_dict.log",
    params:
        extra="",  # optional: extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=10240,
    wrapper:
        "v3.10.2/bio/picard/createsequencedictionary"