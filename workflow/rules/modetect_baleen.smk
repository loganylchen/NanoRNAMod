rule baleen_dataprep:
    input:
        eventalign="results/eventalign/{sample}_baleen.tsv.bz2",
        completion="results/eventalign/{sample}_baleen.tsv.completed",
    output:
        data="results/dataprep/{sample}_baleen_dataprep/eventalign.index",
    params:
        label="{sample}",
    container:
        "docker://btrspg/baleen:dev"
    benchmark:
        "benchmarks/{sample}.baleen_dataprep.txt",
    threads: config['threads']['baleen']
    log:
        out="logs/baleen_dataprep/{sample}.log",
        err="logs/baleen_dataprep/{sample}.error"
    script:
        "../scripts/baleen_dataprep.py"



rule baleen_test:
    input:
        native_eventalign="results/eventalign/{native}_baleen.tsv.bz2",
        native_eventalign_index="results/dataprep/{native}_baleen_dataprep/eventalign.index",
        control_eventalign="results/eventalign/{control}_baleen.tsv.bz2",
        control_eventalign_index="results/dataprep/{control}_baleen_dataprep/eventalign.index",
    output:
        outdir=directory('results/baleen/{native}_{control}/'),
        result='results/baleen/{native}_{control}/transcripts.csv'
    params:
        bedfile=config['target_region']
    container:
        "docker://btrspg/baleen:dev"
    benchmark:
        "benchmarks/{native}_{control}.baleen_test.txt",
    threads: config['threads']['baleen']
    log:
        out="logs/baleen/N_{native}_C_{control}.log",
        err="logs/baleen/N_{native}_C_{control}.error"
    script:
        "../scripts/baleen_mod.py"


