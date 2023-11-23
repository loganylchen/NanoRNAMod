rule tombo_final_lsc:
    input:
        "results/tombo_mod_lsc/{native}_{control}.lsc.tombo.stats"
    output:
        "results/modifications/{native}_{control}/tombo_lsc.tsv.gz"
    log:
        "logs/tombo_lsc_format/{native}_{control}.log"
    shell:
        "touch {output} 2>{log}"

rule tombo_final_msc:
    input:
        "results/tombo_mod_msc/{native}_{control}.msc.tombo.stats"
    output:
        "results/modifications/{native}_{control}/tombo_msc.tsv.gz"
    log:
        "logs/tombo_msc_format/{native}_{control}.log"
    shell:
        "touch {output} 2>{log}"

rule xpore_final:
    input:
        "results/xpore/{native}_{control}/majority_direction_kmer_diffmod.table"
    output:
        "results/modifications/{native}_{control}/xpore.tsv.gz"
    log:
        "logs/xpore_format/{native}_{control}.log"
    shell:
        "touch {output} 2>{log}"

rule nanocompore_final:
    input:
        "results/nanocompore/{native}_{control}/"
    output:
        "results/modifications/{native}_{control}/nanocompore.tsv.gz"
    log:
        "logs/nanocompore_format/{native}_{control}.log"
    shell:
        "touch {output} 2>{log}"

rule m6anet_final:
    input:
        "results/m6anet/{native}/data.site_proba.csv"
    output:
        "results/modifications/{native}_{control}/m6anet.tsv.gz"
    log:
        "logs/m6anet_format/{native}_{control}.log"
    shell:
        "touch {output} 2>{log}"

rule baleen_final:
    input:
        "results/baleen/{native}_{control}/baleen_rmd.done"
    output:
        "results/modifications/{native}_{control}/baleen.tsv.gz"
    log:
        "logs/baleen_format/{native}_{control}.log"
    shell:
        "touch {output} 2>{log}"