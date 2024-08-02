rule post_xpore:
    input:
        "results/xpore/{native}_{control}/diffmod.table",
    output:
        "results/modifications/xpore/{native}_{control}/xpore_results.tsv"
    params:
        tool="xpore",
    log:
        "logs/post_xpore_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_xpore_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule xpore_annotate:
    input:
        transcripts="results/modifications/xpore/{native}_{control}/xpore_results.tsv"
    output:
        result="results/modifications/xpore/{native}_{control}/xpore_annotated_results.tsv"
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config['reference']['transcriptome_gtf'],
    benchmark:
        "benchmarks/{native}_{control}.xpore_annotate.txt",
    threads: 1
    log:
        out="logs/xpore_annotate/N_{native}_C_{control}.log",
        err="logs/xpore_annotate/N_{native}_C_{control}.error"
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column id "
        "--loc-column position "
        "--gtf {param.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_nanocompore:
    input:
        "results/nanocompore/{native}_{control}/nanocompore_results.tsv",
    output:
        "results/modifications/nanocompore/{native}_{control}/nanocompore_results.tsv"
    params:
        tool="nanocompore",
    log:
        "logs/post_nanocompore_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_nanocompore_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule nanocompore_annotate:
    input:
        transcripts="results/modifications/nanocompore/{native}_{control}/nanocompore_results.tsv"
    output:
        result="results/modifications/nanocompore/{native}_{control}/nanocompore_annotated_results.tsv"
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config['reference']['transcriptome_gtf'],
    benchmark:
        "benchmarks/{native}_{control}.nanocompore_annotate.txt",
    threads: 1
    log:
        out="logs/nanocompore_annotate/N_{native}_C_{control}.log",
        err="logs/nanocompore_annotate/N_{native}_C_{control}.error"
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column ref_id "
        "--loc-column pos "
        "--gtf {param.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"



rule post_baleen:
    input:
        'results/baleen/{native}_{control}/transcripts.csv'
    output:
        "results/modifications/baleen/{native}_{control}/baleen_results.tsv"
    params:
        tool="baleen",
    log:
        "logs/post_baleen_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_baleen_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule baleen_annotate:
    input:
        transcripts="results/modifications/baleen/{native}_{control}/baleen_results.tsv"
    output:
        result="results/modifications/baleen/{native}_{control}/baleen_annotated_results.tsv"
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config['reference']['transcriptome_gtf'],
    benchmark:
        "benchmarks/{native}_{control}.baleen_annotate.txt",
    threads: 1
    log:
        out="logs/baleen_annotate/N_{native}_C_{control}.log",
        err="logs/baleen_annotate/N_{native}_C_{control}.error"
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript "
        "--loc-column loc "
        "--gtf {param.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_differr:
    input:
        'results/differr/{native}_{control}/differr.bed'
    output:
        "results/modifications/differr/{native}_{control}/differr_results.tsv"
    params:
        tool="differr",
    log:
        "logs/post_differr_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_differr_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule differr_annotate:
    input:
        transcripts="results/modifications/differr/{native}_{control}/differr_results.tsv"
    output:
        result="results/modifications/differr/{native}_{control}/differr_annotated_results.tsv"
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config['reference']['transcriptome_gtf'],
    benchmark:
        "benchmarks/{native}_{control}.differr_annotate.txt",
    threads: 1
    log:
        out="logs/differr_annotate/N_{native}_C_{control}.log",
        err="logs/differr_annotate/N_{native}_C_{control}.error"
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column chrom "
        "--loc-column start "
        "--gtf {param.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"



rule post_epinano:
    input:
        'results/epinano/{native}_{control}/epinano.delta-sum_err.prediction.csv'
    output:
        "results/modifications/epinano/{native}_{control}/epinano_results.tsv"
    params:
        tool="epinano",
    log:
        "logs/post_epinano_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_epinano_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule epinano_annotate:
    input:
        transcripts="results/modifications/epinano/{native}_{control}/epinano_results.tsv"
    output:
        result="results/modifications/epinano/{native}_{control}/epinano_annotated_results.tsv"
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config['reference']['transcriptome_gtf'],
    benchmark:
        "benchmarks/{native}_{control}.epinano_annotate.txt",
    threads: 1
    log:
        out="logs/epinano_annotate/N_{native}_C_{control}.log",
        err="logs/epinano_annotate/N_{native}_C_{control}.error"
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column chrom "
        "--loc-column pos "
        "--gtf {param.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_eligos2:
    input:
        "results/eligos2/{native}_{control}/{native}_filtered_vs_{control}_filtered_on_{native}_{control}_baseExt0.txt",
    output:
        "results/modifications/eligos2/{native}_{control}/eligos2_results.tsv"
    params:
        tool="eligos2",
    log:
        "logs/post_eligos2_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_eligos2_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule eligos2_annotate:
    input:
        transcripts="results/modifications/eligos2/{native}_{control}/eligos2_results.tsv"
    output:
        result="results/modifications/eligos2/{native}_{control}/eligos2_annotated_results.tsv"
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config['reference']['transcriptome_gtf'],
    benchmark:
        "benchmarks/{native}_{control}.eligos2_annotate.txt",
    threads: 1
    log:
        out="logs/eligos2_annotate/N_{native}_C_{control}.log",
        err="logs/eligos2_annotate/N_{native}_C_{control}.error"
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column chrom "
        "--loc-column start_loc "
        "--gtf {param.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"


rule post_drummer:
    input:
        "results/drummer/{native}_{control}/"
    output:
        "results/modifications/drummer/{native}_{control}/drummer_results.tsv"
    params:
        tool="drummer",
    log:
        "logs/post_drummer_sampled_format/{native}_{control}.log"
    benchmark:
        "benchmarks/{native}_{control}.post_drummer_sampled_format.benchmark.txt"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/format.py"

rule drummer_annotate:
    input:
        transcripts="results/modifications/drummer/{native}_{control}/drummer_results.tsv"
    output:
        result="results/modifications/drummer/{native}_{control}/drummer_annotated_results.tsv"
    container:
        "docker://btrspg/baleen:clean"
    params:
        gtf=config['reference']['transcriptome_gtf'],
    benchmark:
        "benchmarks/{native}_{control}.drummer_annotate.txt",
    threads: 1
    log:
        out="logs/drummer_annotate/N_{native}_C_{control}.log",
        err="logs/drummer_annotate/N_{native}_C_{control}.error"
    shell:
        "Baleen.py annotate "
        "--transcript-mod-file {input.transcripts} "
        "--sep tab "
        "--transcript-column transcript_id "
        "--loc-column transcript_pos "
        "--gtf {param.gtf} "
        "--output {output.result} 2> {log.err} 1> {log.out}"