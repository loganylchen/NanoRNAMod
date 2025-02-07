import glob
import pandas as pd
import sys
import os
from snakemake.utils import validate
from snakemake.logging import logger

validate(config, schema="../schemas/config.schema.yaml")


PROJECT = config["project"]

RESULT_ROOT = f"{PROJECT}/results"
ONLY_FINAL_RESULTS = config["only_final_results"]


def KEEP_OR_NOT(output_file, delete=ONLY_FINAL_RESULTS):
    if delete:
        return temp(output_file)
    else:
        return output_file


# loading samples
samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"SampleName": str}, comment="#")
    .set_index("SampleName", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")

if len(samples["Condition"].unique()) != 2:
    logger.error("Samples should at least come from 2 conditions")
    sys.exit(1)
else:
    native_samples = list(samples.loc[samples["Condition"] == "Native", :].index)
    control_samples = list(samples.loc[samples["Condition"] == "Control", :].index)
    comparisons = [f"{ns}_{cs}" for ns in native_samples for cs in control_samples]


def get_raw_fastq(wildcards):
    if os.path.exists(f"data/{wildcards.sample}/fastq/pass.fq.gz"):
        return f"data/{wildcards.sample}/fastq/pass.fq.gz"
    else:
        return f"data/{wildcards.sample}/fastq/pass.fastq.gz"


def get_raw_blow5(wildcards):
    if os.path.exists(f"data/{wildcards.sample}/blow5/nanopore.blow5"):
        return f"data/{wildcards.sample}/blow5/nanopore.blow5"
    else:
        return f"data/{wildcards.sample}/blow5/nanopore.drs.blow5"


def get_nanocompore_list(sample_list):
    nanocompore_list = [
        f"{RESULT_ROOT}/nanocompore_eventalign_collapse/{sample}/{sample}_eventalign_collapse.tsv"
        for sample in sample_list
    ]
    return ",".join(nanocompore_list)


def get_final_output():
    tools = [tool for tool in config["tools"] if config["tools"][tool]["activate"]]
    final_output = [f"{RESULT_ROOT}/depth_table.tsv"]
    final_output += expand(
        "{RESULT_ROOT}/quantification/{sample}.tx_counts.tsv",
        sample=list(samples.index),
        RESULT_ROOT=RESULT_ROOT,
    )
    final_output += expand(
        "{RESULT_ROOT}/polya/{sample}.genome.nanopolish.polya.tsv.gz",
        sample=list(samples.index),
        RESULT_ROOT=RESULT_ROOT,
    )
    final_output += expand(
        "{RESULT_ROOT}/polya/{sample}.transcriptome.nanopolish.polya.tsv.gz",
        sample=list(samples.index),
        RESULT_ROOT=RESULT_ROOT,
    )
    if config["qc"]:
        final_output += expand(
            "{RESULT_ROOT}/qc/{sample}/{sample}_rnaseq.pdf",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        final_output += expand(
            "{RESULT_ROOT}/qc/{sample}/NanoStats.txt",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        final_output += expand(
            "{RESULT_ROOT}/qc/{sample}/{sample}_samtools_genome_stats.txt",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        final_output += expand(
            "{RESULT_ROOT}/qc/{sample}/{sample}_samtools_transcriptome_stats.txt",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        final_output += expand(
            "{RESULT_ROOT}/variants/{sample}.vcf",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )

    # For some small dataset (on limited transcripts), sampling may be a good choice, while for other big datasets, it may not be necessary

    if "baleen" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/baleen/{comp}/baleen_results.tsv",
            comp=comparisons,
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/baleen/{comp}/baleen_annotated_results.tsv",
                comp=comparisons,
                RESULT_ROOT=RESULT_ROOT,
            )
    if "nanocompore" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/nanocompore/{comp}/nanocompore_results.tsv",
            comp=comparisons,
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/nanocompore/{comp}/nanocompore_annotated_results.tsv",
                comp=comparisons,
                RESULT_ROOT=RESULT_ROOT,
            )
    if "xpore" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/xpore/{comp}/xpore_results.tsv",
            comp=comparisons,
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/xpore/{comp}/xpore_annotated_results.tsv",
                RESULT_ROOT=RESULT_ROOT,
                comp=comparisons,
            )
    if "differr" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/differr/{comp}/differr_results.tsv",
            comp=comparisons,
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/differr/{comp}/differr_annotated_results.tsv",
                comp=comparisons,
                RESULT_ROOT=RESULT_ROOT,
            )
    if "drummer" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/drummer/{comp}/drummer_results.tsv",
            comp=comparisons,
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/drummer/{comp}/drummer_annotated_results.tsv",
                comp=comparisons,
                RESULT_ROOT=RESULT_ROOT,
            )
    if "eligos2" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/eligos2/{comp}/eligos2_results.tsv",
            comp=comparisons,
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/eligos2/{comp}/eligos2_annotated_results.tsv",
                comp=comparisons,
                RESULT_ROOT=RESULT_ROOT,
            )
    if "epinano" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/epinano/{comp}/epinano_results.tsv",
            comp=comparisons,
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/epinano/{comp}/epinano_annotated_results.tsv",
                comp=comparisons,
                RESULT_ROOT=RESULT_ROOT,
            )
    return final_output
