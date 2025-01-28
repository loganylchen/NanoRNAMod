rule post_xpore:
    input:
        "{project}/results/xpore/{native}_{control}/diffmod.table",
    output:
        "{project}/results/modifications/xpore/{native}_{control}/xpore_results.tsv",
    params:
        tool="xpore",
    log:
        "logs/{project}/post_xpore_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_xpore_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: 1
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule xpore_annotate:
    input:
        transcripts="{project}/results/modifications/xpore/{native}_{control}/xpore_results.tsv",
    output:
        result="{project}/results/modifications/xpore/{native}_{control}/xpore_annotated_results.tsv",
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{native}_{control}.xpore_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: 1
    log:
        out="logs/{project}/xpore_annotate/N_{native}_C_{control}.log",
        err="logs/{project}/xpore_annotate/N_{native}_C_{control}.error",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column id "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_nanocompore:
    input:
        "{project}/results/nanocompore/{native}_{control}/nanocompore_results.tsv",
        "{project}/results/nanocompore/{native}_{control}",
    output:
        "{project}/results/modifications/nanocompore/{native}_{control}/nanocompore_results.tsv",
    params:
        tool="nanocompore",
    priority: 20
    log:
        "logs/{project}/post_nanocompore_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_nanocompore_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    threads: 1
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule nanocompore_annotate:
    input:
        transcripts="{project}/results/modifications/nanocompore/{native}_{control}/nanocompore_results.tsv",
    output:
        result="{project}/results/modifications/nanocompore/{native}_{control}/nanocompore_annotated_results.tsv",
    container:
        "docker://btrspg/baleen:clean"
    priority: 20
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.nanocompore_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: 1
    log:
        out="logs/{project}/nanocompore_annotate/N_{native}_C_{control}.log",
        err="logs/{project}/nanocompore_annotate/N_{native}_C_{control}.error",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column ref_id "
        "--loc-column pos "
        "--gtf {params.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_baleen:
    input:
        "{project}/results/baleen/{native}_{control}/transcript_mod_results.csv",
    output:
        "{project}/results/modifications/baleen/{native}_{control}/baleen_results.tsv",
    params:
        tool="baleen",
    priority: 20
    log:
        "logs/{project}/post_baleen_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_baleen_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    threads: 1
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule baleen_annotate:
    input:
        transcripts="{project}/results/modifications/baleen/{native}_{control}/baleen_results.tsv",
    output:
        result="{project}/results/modifications/baleen/{native}_{control}/baleen_annotated_results.tsv",
    container:
        "docker://btrspg/baleen:clean"
    priority: 20
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.baleen_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: 1
    log:
        out="logs/{project}/baleen_annotate/N_{native}_C_{control}.log",
        err="logs/{project}/baleen_annotate/N_{native}_C_{control}.error",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_differr:
    input:
        "{project}/results/differr/{native}_{control}/differr.bed",
    output:
        "{project}/results/modifications/differr/{native}_{control}/differr_results.tsv",
    params:
        tool="differr",
    priority: 20
    log:
        "logs/{project}/post_differr_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_differr_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    threads: 1
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule differr_annotate:
    input:
        transcripts="{project}/results/modifications/differr/{native}_{control}/differr_results.tsv",
    output:
        result="{project}/results/modifications/differr/{native}_{control}/differr_annotated_results.tsv",
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{native}_{control}.differr_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: 1
    log:
        out="logs/{project}/differr_annotate/N_{native}_C_{control}.log",
        err="logs/{project}/differr_annotate/N_{native}_C_{control}.error",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column chrom "
        "--loc-column start "
        "--gtf {params.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_epinano:
    input:
        "{project}/results/epinano/{native}_{control}/epinano.delta-sum_err.prediction.csv",
        "{project}/results/epinano/{native}_{control}",
    output:
        "{project}/results/modifications/epinano/{native}_{control}/epinano_results.tsv",
    params:
        tool="epinano",
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: 1
    log:
        "logs/{project}/post_epinano_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_epinano_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule epinano_annotate:
    input:
        transcripts="{project}/results/modifications/epinano/{native}_{control}/epinano_results.tsv",
    output:
        result="{project}/results/modifications/epinano/{native}_{control}/epinano_annotated_results.tsv",
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.epinano_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: 1
    log:
        out="logs/{project}/epinano_annotate/N_{native}_C_{control}.log",
        err="logs/{project}/epinano_annotate/N_{native}_C_{control}.error",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column chrom "
        "--loc-column pos "
        "--gtf {params.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_eligos2:
    input:
        "{project}/results/eligos2/{native}_{control}/{native}_filtered_vs_{control}_filtered_on_{native}_{control}_baseExt0.txt",
        "{project}/results/eligos2/{native}_{control}",
    output:
        "{project}/results/modifications/eligos2/{native}_{control}/eligos2_results.tsv",
    params:
        tool="eligos2",
    resources:
        mem_mb=1024 * 50,
    threads: 1
    priority: 20
    log:
        "logs/{project}/post_eligos2_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_eligos2_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule eligos2_annotate:
    input:
        transcripts="{project}/results/modifications/eligos2/{native}_{control}/eligos2_results.tsv",
    output:
        result="{project}/results/modifications/eligos2/{native}_{control}/eligos2_annotated_results.tsv",
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.eligos2_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: 1
    priority: 20
    log:
        out="logs/{project}/eligos2_annotate/N_{native}_C_{control}.log",
        err="logs/{project}/eligos2_annotate/N_{native}_C_{control}.error",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column chrom "
        "--loc-column start_loc "
        "--gtf {params.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_drummer:
    input:
        "{project}/results/drummer/{native}_{control}/",
    output:
        "{project}/results/modifications/drummer/{native}_{control}/drummer_results.tsv",
    params:
        tool="drummer",
    resources:
        mem_mb=1024 * 50,
    threads: 1
    priority: 20
    log:
        "logs/{project}/post_drummer_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_drummer_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"


rule drummer_annotate:
    input:
        transcripts="{project}/results/modifications/drummer/{native}_{control}/drummer_results.tsv",
    output:
        result="{project}/results/modifications/drummer/{native}_{control}/drummer_annotated_results.tsv",
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.drummer_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: 1
    priority: 20
    log:
        out="logs/{project}/drummer_annotate/N_{native}_C_{control}.log",
        err="logs/{project}/drummer_annotate/N_{native}_C_{control}.error",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript_id "
        "--loc-column transcript_pos "
        "--gtf {params.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"
