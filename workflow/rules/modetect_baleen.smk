rule baleen_dataprep:
    input:
        eventalign="results/eventalign/{sample}_full.tsv.bz2",
        completion="results/eventalign/{sample}_full.tsv.completed",
    output:
        data="results/dataprep/{sample}_baleen_dataprep/eventalign.index",
        dir=directory("results/dataprep/{sample}_baleen_dataprep"),
    params:
        label="{sample}",
        use_mem=config['baleen']['use_mem'],
    container:
        "docker://btrspg/baleen:dev"
    benchmark:
        "benchmarks/{sample}.baleen_dataprep.txt",
    threads: config['threads']['baleen']
    log:
        out="logs/baleen_dataprep/{sample}.log",
        err="logs/baleen_dataprep/{sample}.error"
    shell:
        "Baleen.py dataprep "
        "--eventalign {input.eventalign} "
        "--output-dir {output.dir} "
        "--threads {threads} 1> {log.out} 2> {log.err}"


rule baleen_dataprep_sampled:
    input:
        eventalign="results/eventalign/{sample}_full_{sample_size}_{n}.tsv.bz2",
        completion="results/eventalign/{sample}_full_{sample_size}_{n}.tsv.completed",
    output:
        data="results/dataprep/{sample}_baleen_dataprep_{sample_size}_{n}/eventalign.index",
        dir=directory("results/dataprep/{sample}_baleen_dataprep_{sample_size}_{n}"),
    params:
        label="{sample}",
        use_mem=config['baleen']['use_mem'],
    container:
        "docker://btrspg/baleen:dev"
    benchmark:
        "benchmarks/{sample}.baleen_dataprep_{sample_size}_{n}.txt",
    threads: config['threads']['baleen']
    log:
        out="logs/baleen_dataprep/{sample}_{sample_size}_{n}.log",
        err="logs/baleen_dataprep/{sample}_{sample_size}_{n}.error"
    shell:
        "Baleen.py dataprep "
        "--eventalign {input.eventalign} "
        "--output-dir {output.dir} "
        "--threads {threads} 1> {log.out} 2> {log.err}"


rule baleen_test:
    input:
        native_dataprep="results/dataprep/{native}_baleen_dataprep/",
        native_eventalign_index="results/dataprep/{native}_baleen_dataprep/eventalign.index",
        control_dataprep="results/dataprep/{control}_baleen_dataprep/",
        control_eventalign_index="results/dataprep/{control}_baleen_dataprep/eventalign.index",
    output:
        result='results/baleen/{native}_{control}/transcripts.csv'
    params:
        bedfile=config['target_region'],
        params=config['params']['baleen_modcall'],
        dir="results/baleen/{native}_{control}/",
    container:
        "docker://btrspg/baleen:dev"
    benchmark:
        "benchmarks/{native}_{control}.baleen_test.txt",
    threads: config['threads']['baleen']
    log:
        out="logs/baleen/N_{native}_C_{control}.log",
        err="logs/baleen/N_{native}_C_{control}.error"
    shell:
        "Baleen.py modcall "
        "--native-dataprep {input.native_dataprep} "
        "--control-dataprep {input.control_dataprep} "
        "{params.params} "
        "--output-dir {params.dir} 1>{log.out} 2> {log.err}"


rule baleen_test_sampled:
    input:
        native_dataprep="results/dataprep/{native}_baleen_dataprep_{sample_size}_{n}",
        native_eventalign_index="results/dataprep/{native}_baleen_dataprep_{sample_size}_{n}/eventalign.index",
        control_dataprep="results/dataprep/{control}_baleen_dataprep_{sample_size}_{n}",
        control_eventalign_index="results/dataprep/{control}_baleen_dataprep_{sample_size}_{n}/eventalign.index",
    output:
        result='results/baleen/{native}_{control}-{sample_size}-{n}/transcripts.csv'
    params:
        bedfile=config['target_region'],
        params=config['params']['baleen_modcall'],
        dir= "results/baleen/{native}_{control}-{sample_size}-{n}/",
    container:
        "docker://btrspg/baleen:dev"
    benchmark:
        "benchmarks/{native}_{control}_{sample_size}_{n}.baleen_test.txt",
    threads: config['threads']['baleen']
    log:
        out="logs/baleen/N_{native}_C_{control}_{sample_size}_{n}.log",
        err="logs/baleen/N_{native}_C_{control}_{sample_size}_{n}.error"
    shell:
        "Baleen.py modcall "
        "--native-dataprep {input.native_dataprep} "
        "--control-dataprep {input.control_dataprep} "
        "{params.params} "
        "--output-dir {params.dir} 1>{log.out} 2>{log.err} "