rule post_tandemmod:
    input:
        predictions="{project}/results/dataprep/{sample}_tandemmod_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_tandemmod_dataprep"
    output:
        "{project}/results/modifications/tandemmod/{sample}/tandemmod_results.tsv",
    params:
        tool="tandemmod",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_tandemmod_sampled_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_tandemmod_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule tandemmod_annotate:
    input:
        transcripts="{project}/results/modifications/tandemmod/{sample}/tandemmod_results.tsv",
    output:
        result="{project}/results/modifications/tandemmod/{sample}/tandemmod_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{sample}.tandemmod_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/tandemmod_annotate/{sample}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "> {log} 2>&1"


rule post_directrm:
    input:
        predictions="{project}/results/dataprep/{sample}_directrm_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_directrm_dataprep"
    output:
        "{project}/results/modifications/directrm/{sample}/directrm_results.tsv",
    params:
        tool="directrm",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_directrm_sampled_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_directrm_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule directrm_annotate:
    input:
        transcripts="{project}/results/modifications/directrm/{sample}/directrm_results.tsv",
    output:
        result="{project}/results/modifications/directrm/{sample}/directrm_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{sample}.directrm_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/directrm_annotate/{sample}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "> {log} 2>&1"


rule post_m6atm:
    input:
        predictions="{project}/results/dataprep/{sample}_m6atm_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_m6atm_dataprep"
    output:
        "{project}/results/modifications/m6atm/{sample}/m6atm_results.tsv",
    params:
        tool="m6atm",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_m6atm_sampled_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_m6atm_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule m6atm_annotate:
    input:
        transcripts="{project}/results/modifications/m6atm/{sample}/m6atm_results.tsv",
    output:
        result="{project}/results/modifications/m6atm/{sample}/m6atm_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{sample}.m6atm_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/m6atm_annotate/{sample}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "> {log} 2>&1"


rule post_rnano:
    input:
        predictions="{project}/results/dataprep/{sample}_rnano_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_rnano_dataprep"
    output:
        "{project}/results/modifications/rnano/{sample}/rnano_results.tsv",
    params:
        tool="rnano",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_rnano_sampled_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_rnano_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule rnano_annotate:
    input:
        transcripts="{project}/results/modifications/rnano/{sample}/rnano_results.tsv",
    output:
        result="{project}/results/modifications/rnano/{sample}/rnano_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{sample}.rnano_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/rnano_annotate/{sample}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "> {log} 2>&1"


rule post_xpore:
    input:
        "{project}/results/xpore/{native}_{control}/diffmod.table",
    output:
        "{project}/results/modifications/xpore/{native}_{control}/xpore_results.tsv",
    params:
        tool="xpore",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_xpore_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_xpore_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule xpore_annotate:
    input:
        transcripts="{project}/results/modifications/xpore/{native}_{control}/xpore_results.tsv",
    output:
        result="{project}/results/modifications/xpore/{native}_{control}/xpore_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{native}_{control}.xpore_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/xpore_annotate/N_{native}_C_{control}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column id "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} > {log} 2>&1"


rule post_nanocompore:
    input:
        "{project}/results/nanocompore/{native}_{control}/nanocompore_results.tsv",
        "{project}/results/nanocompore/{native}_{control}",
    output:
        "{project}/results/modifications/nanocompore/{native}_{control}/nanocompore_results.tsv",
    params:
        tool="nanocompore",
    container:
        get_container("python3")
    priority: 20
    log:
        "logs/{project}/post_nanocompore_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_nanocompore_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule nanocompore_annotate:
    input:
        transcripts="{project}/results/modifications/nanocompore/{native}_{control}/nanocompore_results.tsv",
    output:
        result="{project}/results/modifications/nanocompore/{native}_{control}/nanocompore_annotated_results.tsv",
    container:
        get_container("baleen")
    priority: 20
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.nanocompore_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/nanocompore_annotate/N_{native}_C_{control}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column ref_id "
        "--loc-column pos "
        "--gtf {params.gtf} "
        "--output {output.result} > {log} 2>&1"


rule post_baleen:
    input:
        "{project}/results/baleen/{native}_{control}/transcript_mod_results.csv",
    output:
        "{project}/results/modifications/baleen/{native}_{control}/baleen_results.tsv",
    params:
        tool="baleen",
    container:
        get_container("python3")
    priority: 20
    log:
        "logs/{project}/post_baleen_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_baleen_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule baleen_annotate:
    input:
        transcripts="{project}/results/modifications/baleen/{native}_{control}/baleen_results.tsv",
    output:
        result="{project}/results/modifications/baleen/{native}_{control}/baleen_annotated_results.tsv",
    container:
        get_container("baleen")
    priority: 20
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.baleen_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/baleen_annotate/N_{native}_C_{control}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} > {log} 2>&1"


rule post_differr:
    input:
        "{project}/results/differr/{native}_{control}/differr.bed",
    output:
        "{project}/results/modifications/differr/{native}_{control}/differr_results.tsv",
    params:
        tool="differr",
    container:
        get_container("python3")
    priority: 20
    log:
        "logs/{project}/post_differr_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_differr_sampled_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule differr_annotate:
    input:
        transcripts="{project}/results/modifications/differr/{native}_{control}/differr_results.tsv",
    output:
        result="{project}/results/modifications/differr/{native}_{control}/differr_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{native}_{control}.differr_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/differr_annotate/N_{native}_C_{control}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column chrom "
        "--loc-column start "
        "--gtf {params.gtf} "
        "--output {output.result} > {log} 2>&1"


rule post_epinano:
    input:
        "{project}/results/epinano/{native}_{control}/epinano.delta-sum_err.prediction.csv",
        "{project}/results/epinano/{native}_{control}",
    output:
        "{project}/results/modifications/epinano/{native}_{control}/epinano_results.tsv",
    params:
        tool="epinano",
    container:
        get_container("python3")
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    log:
        "logs/{project}/post_epinano_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_epinano_sampled_format.benchmark.txt"
    script:
        "../scripts/format.py"


rule epinano_annotate:
    input:
        transcripts="{project}/results/modifications/epinano/{native}_{control}/epinano_results.tsv",
    output:
        result="{project}/results/modifications/epinano/{native}_{control}/epinano_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.epinano_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/epinano_annotate/N_{native}_C_{control}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column chrom "
        "--loc-column pos "
        "--gtf {params.gtf} "
        "--output {output.result} > {log} 2>&1"


rule post_eligos2:
    input:
        "{project}/results/eligos2/{native}_{control}/{native}_vs_{control}_on_{native}_{control}_baseExt0.txt",
        "{project}/results/eligos2/{native}_{control}",
    output:
        "{project}/results/modifications/eligos2/{native}_{control}/eligos2_results.tsv",
    params:
        tool="eligos2",
    container:
        get_container("python3")
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("python3", 1)
    priority: 20
    log:
        "logs/{project}/post_eligos2_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_eligos2_sampled_format.benchmark.txt"
    script:
        "../scripts/format.py"


rule eligos2_annotate:
    input:
        transcripts="{project}/results/modifications/eligos2/{native}_{control}/eligos2_results.tsv",
    output:
        result="{project}/results/modifications/eligos2/{native}_{control}/eligos2_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.eligos2_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    priority: 20
    log:
        "logs/{project}/eligos2_annotate/N_{native}_C_{control}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column chrom "
        "--loc-column start_loc "
        "--gtf {params.gtf} "
        "--output {output.result} > {log} 2>&1"


rule post_drummer:
    input:
        "{project}/results/drummer/{native}_{control}/",
    output:
        "{project}/results/modifications/drummer/{native}_{control}/drummer_results.tsv",
    params:
        tool="drummer",
    container:
        get_container("python3")
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("python3", 1)
    priority: 20
    log:
        "logs/{project}/post_drummer_sampled_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_drummer_sampled_format.benchmark.txt"
    script:
        "../scripts/format.py"


rule drummer_annotate:
    input:
        transcripts="{project}/results/modifications/drummer/{native}_{control}/drummer_results.tsv",
    output:
        result="{project}/results/modifications/drummer/{native}_{control}/drummer_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    benchmark:
        "benchmarks/{project}/{native}_{control}.drummer_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    priority: 20
    log:
        "logs/{project}/drummer_annotate/N_{native}_C_{control}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript_id "
        "--loc-column transcript_pos "
        "--gtf {params.gtf} "
        "--output {output.result} > {log} 2>&1"


rule post_psipore:
    input:
        "{project}/results/psipore/{native}_{control}/psipore_results.tsv",
    output:
        "{project}/results/modifications/psipore/{native}_{control}/psipore_results.tsv",
    params:
        tool="psipore",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_psipore_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_psipore_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule psipore_annotate:
    input:
        transcripts="{project}/results/modifications/psipore/{native}_{control}/psipore_results.tsv",
    output:
        result="{project}/results/modifications/psipore/{native}_{control}/psipore_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{native}_{control}.psipore_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/psipore_annotate/N_{native}_C_{control}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "> {log} 2>&1"


rule post_nanopsu:
    input:
        predictions="{project}/results/dataprep/{sample}_nanopsu_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_nanopsu_dataprep",
    output:
        "{project}/results/modifications/nanopsu/{sample}/nanopsu_results.tsv",
    params:
        tool="nanopsu",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_nanopsu_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_nanopsu_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule nanopsu_annotate:
    input:
        transcripts="{project}/results/modifications/nanopsu/{sample}/nanopsu_results.tsv",
    output:
        result="{project}/results/modifications/nanopsu/{sample}/nanopsu_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{sample}.nanopsu_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/nanopsu_annotate/{sample}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "> {log} 2>&1"


rule post_nanomud:
    input:
        predictions="{project}/results/dataprep/{sample}_nanomud_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_nanomud_dataprep",
    output:
        "{project}/results/modifications/nanomud/{sample}/nanomud_results.tsv",
    params:
        tool="nanomud",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_nanomud_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_nanomud_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/format.py"


rule nanomud_annotate:
    input:
        transcripts="{project}/results/modifications/nanomud/{sample}/nanomud_results.tsv",
    output:
        result="{project}/results/modifications/nanomud/{sample}/nanomud_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{sample}.nanomud_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/nanomud_annotate/{sample}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "> {log} 2>&1"


rule post_penguin:
    input:
        predictions="{project}/results/dataprep/{sample}_penguin_dataprep/predictions.tsv",
        directory="{project}/results/dataprep/{sample}_penguin_dataprep",
    output:
        "{project}/results/modifications/penguin/{sample}/penguin_results.tsv",
    params:
        tool="penguin",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_penguin_format/{sample}.log",
    benchmark:
        "benchmarks/{project}/{sample}.post_penguin_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("default", 1)
    script:
        "../scripts/format.py"


rule penguin_annotate:
    input:
        transcripts="{project}/results/modifications/penguin/{sample}/penguin_results.tsv",
    output:
        result="{project}/results/modifications/penguin/{sample}/penguin_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{sample}.penguin_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/penguin_annotate/{sample}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "> {log} 2>&1"


# pybaleen post-processing
rule post_pybaleen:
    input:
        result="{project}/results/pybaleen/{native}_{control}/site_results.tsv",
    output:
        result="{project}/results/modifications/pybaleen/{native}_{control}/pybaleen_results.tsv",
    container:
        get_container("python3")
    log:
        "logs/{project}/post_pybaleen_format/{native}_{control}.log",
    benchmark:
        "benchmarks/{project}/{native}_{control}.post_pybaleen_format.benchmark.txt"
    resources:
        mem_mb=1024 * 50,
    priority: 20
    threads: get_threads("python3", 1)
    script:
        "../scripts/pybaleen_postprocess.py"


rule pybaleen_annotate:
    input:
        transcripts="{project}/results/modifications/pybaleen/{native}_{control}/pybaleen_results.tsv",
    output:
        result="{project}/results/modifications/pybaleen/{native}_{control}/pybaleen_annotated_results.tsv",
    container:
        get_container("baleen")
    params:
        gtf=config["reference"]["transcriptome_gtf"],
    priority: 20
    benchmark:
        "benchmarks/{project}/{native}_{control}.pybaleen_annotate.txt"
    resources:
        mem_mb=1024 * 50,
    threads: get_threads("baleen", 1)
    log:
        "logs/{project}/pybaleen_annotate/N_{native}_C_{control}.log",
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column position "
        "--gtf {params.gtf} "
        "--output {output.result} "
        "> {log} 2>&1"
