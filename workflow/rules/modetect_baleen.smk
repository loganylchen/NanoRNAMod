rule baleen_dataprep:
    input:
        eventalign="{project}/results/eventalign/{sample}_full.tsv.bz2",
        completion="{project}/results/eventalign/{sample}_full.tsv.completed",
    output:
        data="{project}/results/dataprep/{sample}_baleen_dataprep/eventalign.index",
        dir=KEEP_OR_NOT(
            directory("{project}/results/dataprep/{sample}_baleen_dataprep")
        ),
    container:
        "docker://btrspg/baleen:clean"
    benchmark:
        "benchmarks/{project}/{sample}.baleen_dataprep.txt"
    threads: config["threads"]["baleen"]
    resources:
        mem_mb=1024 * 450,
    log:
        out="logs/{project}/baleen_dataprep/{sample}.log",
        err="logs/{project}/baleen_dataprep/{sample}.error",
    shell:
        "Baleen.py dataprep "
        "--eventalign {input.eventalign} "
        "--output-dir {output.dir} "
        "--threads {threads} 1> {log.out} 2> {log.err}"


rule baleen_modcall:
    input:
        native_dataprep="{project}/results/dataprep/{native}_baleen_dataprep/",
        native_eventalign_index="{project}/results/dataprep/{native}_baleen_dataprep/eventalign.index",
        control_dataprep="{project}/results/dataprep/{control}_baleen_dataprep/",
        control_eventalign_index="{project}/results/dataprep/{control}_baleen_dataprep/eventalign.index",
    output:
        result=KEEP_OR_NOT(
            "{project}/results/baleen/{native}_{control}/transcript_mod_results.csv"
        ),
    params:
        bedfile=config["target_region"],
        params=config["params"]["baleen_modcall"],
        output_dir="{project}/results/baleen/{native}_{control}/",
    container:
        "docker://btrspg/baleen:clean"
    benchmark:
        "benchmarks/{project}/{native}_{control}.baleen_modcall.txt"
    threads: config["threads"]["baleen"]
    resources:
        mem_mb=1024 * 650,
    log:
        out="logs/{project}/baleen_modcall/N_{native}_C_{control}.log",
        err="logs/{project}/baleen_modcall/N_{native}_C_{control}.error",
    shell:
        "Baleen.py modcall "
        "--native-dataprep {input.native_dataprep} "
        "--control-dataprep {input.control_dataprep} "
        "{params.params} "
        "--output-dir {params.output_dir} 1>{log.out} 2> {log.err}"


rule baleen_postcall:
    input:
        modcall="{project}/results/baleen/{native}_{control}/modcall_sm",
    output:
        result="{project}/results/baleen/{native}_{control}/transcripts.csv",
    params:
        params=config["params"]["baleen_postcall"],
        dir="{project}/results/baleen/{native}_{control}/",
    container:
        "docker://btrspg/baleen:clean"
    benchmark:
        "benchmarks/{project}/{native}_{control}.baleen_postcall.txt"
    threads: config["threads"]["baleen"]
    resources:
        mem_mb=1024 * 250,
    log:
        out="logs/{project}/baleen_postcall/N_{native}_C_{control}.log",
        err="logs/{project}/baleen_postcall/N_{native}_C_{control}.error",
    shell:
        "Baleen.py postcall "
        "--modcall-sm-dir {input.modcall} "
        "--threads {threads} "
        "{params.params} "
        "--output-dir {params.dir} 2> {log.err} 1> {log.out}"
