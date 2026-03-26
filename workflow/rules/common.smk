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


PER_COMPARISON_TOOLS = ["xpore", "nanocompore", "baleen", "differr", "drummer",
                        "eligos2", "epinano", "psipore", "pybaleen"]


def get_active_comparison_tools():
    return [t for t in config["tools"]
            if config["tools"][t]["activate"] and t in PER_COMPARISON_TOOLS]


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
    final_output = [f"{RESULT_ROOT}/workflow_version.json"]  # Always record version first
    final_output += [f"{RESULT_ROOT}/depth_table.tsv"]
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
    # Benchmarking outputs - resource metrics (shared across all tools)
    final_output += [f"{RESULT_ROOT}/benchmarks/resource_summary.tsv"]
    final_output += [f"{RESULT_ROOT}/benchmarks/resource_by_tool.tsv"]
    final_output += [f"{RESULT_ROOT}/benchmarks/resource_by_tool.pdf"]
    # Aggregated accuracy benchmarking outputs (produced by accuracy_benchmark.py)
    if config.get("benchmark", {}).get("truth_set", ""):
        active_tools = get_active_comparison_tools()

        # Per-tool coverage
        final_output += expand(
            f"{RESULT_ROOT}/benchmarks/coverage/{{comp}}/{{tool}}_covered.tsv",
            comp=comparisons, tool=active_tools,
        )
        final_output += expand(
            f"{RESULT_ROOT}/benchmarks/coverage/{{comp}}/union.tsv",
            comp=comparisons,
        )

        # Per-tool native and fair evaluation
        final_output += expand(
            f"{RESULT_ROOT}/benchmarks/native/{{tool}}/{{comp}}/best_metrics.tsv",
            tool=active_tools, comp=comparisons,
        )
        final_output += expand(
            f"{RESULT_ROOT}/benchmarks/fair/{{tool}}/{{comp}}/best_metrics.tsv",
            tool=active_tools, comp=comparisons,
        )

        # Aggregated summaries (moved to aggregated/ subdir)
        final_output += [f"{RESULT_ROOT}/benchmarks/aggregated/accuracy_summary.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/aggregated/accuracy_summary_overall.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/aggregated/accuracy_summary_by_comparison.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/aggregated/by_tool.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/aggregated/best_scores.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/aggregated/called_sites_by_comparison.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/aggregated/called_sites_summary.tsv"]
        # Visualization and detailed reports (shared across all tools)
        final_output += [f"{RESULT_ROOT}/benchmarks/viz/benchmark_report.html"]
        final_output += [f"{RESULT_ROOT}/benchmarks/benchmark_report.pdf"]
        final_output += [f"{RESULT_ROOT}/benchmarks/optimal_thresholds.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/detailed_predictions.tsv"]
        # Multi-threshold evaluation outputs
        final_output += [f"{RESULT_ROOT}/benchmarks/threshold_evaluation.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/optimal_thresholds_detailed.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/score_distributions.tsv"]
        # Negative control strategies (k-mer and same-base negatives)
        final_output += [f"{RESULT_ROOT}/benchmarks/accuracy_summary_kmer_negatives.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/kmer_negative_sites.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/accuracy_summary_same_base_negatives.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/same_base_negative_sites.tsv"]
        # Statistical analysis (bootstrap CIs, significance tests, effect sizes)
        final_output += [f"{RESULT_ROOT}/benchmarks/statistics/bootstrap_ci.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/statistics/significance_tests.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/statistics/fdr_corrected.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/statistics/effect_sizes.tsv"]
        # Sensitivity analysis (coverage, score distribution, threshold robustness)
        final_output += [f"{RESULT_ROOT}/benchmarks/sensitivity/coverage_analysis.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/sensitivity/score_distribution.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/sensitivity/threshold_robustness.tsv"]
        # Publication figures (ggplot2, Nature/Cell/Science themes)
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/main/fig1_overall_accuracy.pdf"]
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/main/fig2_pr_curves.pdf"]
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/main/fig3_roc_curves.pdf"]
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/main/fig4_f1_comparison.pdf"]
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/main/fig5_coverage_sensitivity.pdf"]
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/main/fig6_resource_usage.pdf"]
        # Supplementary figures
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/supplementary/sfig1_per_comparison.pdf"]
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/supplementary/sfig2_score_distributions.pdf"]
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/supplementary/sfig3_threshold_robustness.pdf"]
        final_output += [f"{RESULT_ROOT}/benchmarks/figures/supplementary/sfig4_effect_sizes.pdf"]
        # Figure source data (for manuscript submission)
        final_output += [f"{RESULT_ROOT}/benchmarks/data/fig1_source_data.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/data/fig2_source_data.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/data/fig3_source_data.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/data/fig4_source_data.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/data/fig5_source_data.tsv"]
        final_output += [f"{RESULT_ROOT}/benchmarks/data/fig6_source_data.tsv"]
    return final_output


##### version tracking rule #####

rule workflow_version:
    """
    Record workflow version information for reproducibility.

    Creates a JSON file with:
    - workflow version (git commit SHA or snakedeploy version)
    - project name
    - timestamp of workflow run

    Uses WORKFLOW_VERSION from Snakefile (computed at startup).
    """
    output:
        touch(f"{RESULT_ROOT}/workflow_version.json")
    run:
        import json
        from datetime import datetime

        # Use WORKFLOW_VERSION from the Snakefile (passed via globals)
        version = WORKFLOW_VERSION

        version_info = {
            "workflow": "NanoRNAMod",
            "version": version,
            "project": PROJECT,
            "config_file": str(configfile if 'configfile' in dir() else config.get("configfile", "unknown")),
            "timestamp": datetime.now().isoformat(),
        }

        with open(output[0], "w") as f:
            json.dump(version_info, f, indent=2)

        logger.info(f"Workflow version saved: {version}")

