rule nanocompore_collapse:
    """
    Decompress the bz2 eventalign output to a tmpdir-local tempfile,
    run nanocompore's collapse, then auto-cleanup via trap. Replaces
    the previous uncompress_eventalign_full → nanocompore_collapse
    two-rule chain.
    """
    input:
        eventalign="{project}/results/eventalign/{sample}_full.tsv.bz2",
        completion="{project}/results/eventalign/{sample}_full.tsv.completed",
    output:
        output=temp(
            "{project}/results/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv"
        ),
    params:
        prefix="{sample}",
        dir="{project}/results/nanocompore_eventalign_collapse/{sample}",
    log:
        "logs/{project}/nanocompore_collapse/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.nanocompore_collapse.benchmark.txt"
    container:
        get_container("nanocompore")
    resources:
        mem_mb = 1024 * 50
    priority: 84
    threads: get_threads("nanocompore", 4)
    shell:
        "mkdir -p {params.dir}; "
        "TMP={params.dir}/{params.prefix}.eventalign.uncompressed.tsv; "
        "trap 'rm -f $TMP' EXIT; "
        "bzip2 -dc {input.eventalign} > $TMP && "
        "nanocompore eventalign_collapse "
        "-i $TMP "
        "--outpath {params.dir} "
        "--outprefix {params.prefix} "
        "--overwrite "
        "--nthreads {threads} 1>{log} 2>&1"


rule nanocompore:
    input:
        control_file="{project}/results/nanocompore_eventalign_collapse/{control}/{control}_eventalign_collapse.tsv",
        native_file="{project}/results/nanocompore_eventalign_collapse/{native}/{native}_eventalign_collapse.tsv",
        reference=config["reference"]["transcriptome_fasta"],
    output:
        output_dir=temp(
            directory("{project}/results/nanocompore/{native}_{control}")
        ),
        output_file=temp("{project}/results/nanocompore/{native}_{control}/nanocompore_results.tsv"),
    params:
        prefix="{native}_{control}",
        extra=config["params"]["nanocompore"],
    container:
        get_container("nanocompore")
    priority: 85
    resources:
        mem_mb = 1024 * 50
    threads: get_threads("nanocompore", 4)
    log:
        "logs/{project}/nanocompore/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.nanocompore.benchmark.txt"
    shell:
        "nanocompore sampcomp "
        "--file_list1 {input.control_file} "
        "--file_list2 {input.native_file} "
        "--label1 Control "
        "--label2 Native "
        "--nthreads {threads} "
        "{params.extra} "
        "--fasta {input.reference} "
        "--outpath {output.output_dir} "
        "--outprefix '' "
        "--overwrite  1>{log} 2>&1; "
        "status=$?; "
        "if [ $status -eq 0 ]; then "
        "  mkdir -p {output.output_dir}; "
        "  [ -f {output.output_file} ] || touch {output.output_file}; "
        "fi; "
        "exit $status"
