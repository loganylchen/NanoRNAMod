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


def get_container(tool_name):
    if "containers" in config and tool_name in config["containers"]:
        container_image = config["containers"][tool_name]
        if container_image and container_image.strip():
            return container_image
    if "containers" in config and "default" in config["containers"]:
        return config["containers"]["default"]
    return "docker://condaforge/mambaforge:22.11.1-4"


def get_threads(tool_name, default=1):
    if "threads" in config and tool_name in config["threads"]:
        return config["threads"][tool_name]
    if "threads" in config and "default" in config["threads"]:
        return config["threads"]["default"]
    return default


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
            "{RESULT_ROOT}/variants/{sample}.bcf",
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
    if "tandemmod" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/tandemmod/{sample}/tandemmod_results.tsv",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/tandemmod/{sample}/tandemmod_annotated_results.tsv",
                sample=list(samples.index),
                RESULT_ROOT=RESULT_ROOT,
            )
    if "directrm" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/directrm/{sample}/directrm_results.tsv",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/directrm/{sample}/directrm_annotated_results.tsv",
                sample=list(samples.index),
                RESULT_ROOT=RESULT_ROOT,
            )
    if "m6atm" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/m6atm/{sample}/m6atm_results.tsv",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/m6atm/{sample}/m6atm_annotated_results.tsv",
                sample=list(samples.index),
                RESULT_ROOT=RESULT_ROOT,
            )
    if "rnano" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/rnano/{sample}/rnano_results.tsv",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/rnano/{sample}/rnano_annotated_results.tsv",
                sample=list(samples.index),
                RESULT_ROOT=RESULT_ROOT,
            )
    if "psipore" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/psipore/{comp}/psipore_results.tsv",
            comp=comparisons,
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/psipore/{comp}/psipore_annotated_results.tsv",
                comp=comparisons,
                RESULT_ROOT=RESULT_ROOT,
            )
    if "nanopsu" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/nanopsu/{sample}/nanopsu_results.tsv",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/nanopsu/{sample}/nanopsu_annotated_results.tsv",
                sample=list(samples.index),
                RESULT_ROOT=RESULT_ROOT,
            )
    if "nanomud" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/nanomud/{sample}/nanomud_results.tsv",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/nanomud/{sample}/nanomud_annotated_results.tsv",
                sample=list(samples.index),
                RESULT_ROOT=RESULT_ROOT,
            )
    if "penguin" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/penguin/{sample}/penguin_results.tsv",
            sample=list(samples.index),
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/penguin/{sample}/penguin_annotated_results.tsv",
                sample=list(samples.index),
                RESULT_ROOT=RESULT_ROOT,
            )
    if "pybaleen" in tools:
        final_output += expand(
            "{RESULT_ROOT}/modifications/pybaleen/{comp}/pybaleen_results.tsv",
            comp=comparisons,
            RESULT_ROOT=RESULT_ROOT,
        )
        if os.path.exists(config["reference"]["transcriptome_gtf"]):
            final_output += expand(
                "{RESULT_ROOT}/modifications/pybaleen/{comp}/pybaleen_annotated_results.tsv",
                comp=comparisons,
                RESULT_ROOT=RESULT_ROOT,
            )
    # Benchmarking outputs
    final_output += [f"{RESULT_ROOT}/benchmarks/resource_summary.tsv"]
    if config.get("benchmark", {}).get("truth_set", ""):
        final_output += [f"{RESULT_ROOT}/benchmarks/accuracy_summary.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/accuracy_summary_overall.tsv"]
        # Visualization and detailed reports
        final_output += [f"{RESULT_ROOT}/benchmarks/viz/benchmark_report.html"]
        final_output += [f"{RESULT_ROOT}/benchmarks/optimal_thresholds.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/detailed_predictions.tsv"]
        # Multi-threshold evaluation outputs
        final_output += [f"{RESULT_ROOT}/benchmarks/threshold_evaluation.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/optimal_thresholds_detailed.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/score_distributions.tsv"]
    return final_output
