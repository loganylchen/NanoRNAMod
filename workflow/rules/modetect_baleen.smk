rule baleen_dataprep:
    input:
        eventalign="results/eventalign/{sample}_full.tsv.bz2",
        completion="results/eventalign/{sample}_full.tsv.completed",
    output:
        data="results/dataprep/{sample}_baleen_dataprep/eventalign.index",
        dir=directory("results/dataprep/{sample}_baleen_dataprep"),
    container:
        "docker://btrspg/baleen:clean"
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



rule baleen_modcall:
    input:
        native_dataprep="results/dataprep/{native}_baleen_dataprep/",
        native_eventalign_index="results/dataprep/{native}_baleen_dataprep/eventalign.index",
        control_dataprep="results/dataprep/{control}_baleen_dataprep/",
        control_eventalign_index="results/dataprep/{control}_baleen_dataprep/eventalign.index",
    output:
        result=directory('results/baleen/{native}_{control}/modcall_sm'),
    params:
        bedfile=config['target_region'],
        params=config['params']['baleen_modcall'],
        dir="results/baleen_modcall/{native}_{control}/",
    container:
        "docker://btrspg/baleen:clean"
    benchmark:
        "benchmarks/{native}_{control}.baleen_modcall.txt",
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



rule baleen_postcall:
    input:
        modcall="results/baleen/{native}_{control}/modcall_sm",
    output:
        result='results/baleen/{native}_{control}/transcripts.csv'
    params:
        params=config['params']['baleen_postcall'],
        dir="results/baleen/{native}_{control}/",
    container:
        "docker://btrspg/baleen:clean"
    benchmark:
        "benchmarks/{native}_{control}.baleen_postcall.txt",
    threads: config['threads']['baleen']
    log:
        "logs/baleen_postcall/N_{native}_C_{control}.log",
    shell:
        "Baleen.py postcall "
        "--modcall-sm-dir {input.modcall} "
        "--threads {threads} "
        "{params.params} "
        "--output-dir {params.dir} 2> {log}"
