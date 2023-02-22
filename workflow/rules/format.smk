rule tombo_final_lsc:
    input:
        "results/tombo_mod_lsc/{native}_{control}.lsc.tombo.stats"
    output:
        "results/modifications/{native}_{control}/tombo_lsc.tsv.gz"
    shell:
        "touch {output}"

rule tombo_final_msc:
    input:
        "results/tombo_mod_msc/{native}_{control}.msc.tombo.stats"
    output:
        "results/modifications/{native}_{control}/tombo_msc.tsv.gz"
    shell:
        "touch {output}"

rule xpore_final:
    input:
        "results/xpore/{native}_{control}/majority_direction_kmer_diffmod.table"
    output:
        "results/modifications/{native}_{control}/xpore.tsv.gz"
    shell:
        "touch {outpu}"

rule nanocompore_final:
    input:
        "results/nanocompore/{native}_{control}/"
    output:
        "results/modifications/{native}_{control}/nanocompore.tsv.gz"
    shell:
        "touch {output}"

rule m6anet_final:
    input:
        "results/m6anet/{native}/data.result.csv.gz"
    output:
        "results/modifications/{native}_{control}/m6anet.tsv.gz"
    shell:
        "touch {output}"