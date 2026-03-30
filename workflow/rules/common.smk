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

PER_SAMPLE_TOOLS = ["tandemmod", "directrm", "m6atm", "rnano", "nanopsu", "nanomud", "penguin"]


def get_active_comparison_tools():
    return [t for t in config["tools"]
            if config["tools"][t]["activate"] and t in PER_COMPARISON_TOOLS]


def get_active_per_sample_tools():
    return [t for t in config["tools"]
            if config["tools"][t]["activate"] and t in PER_SAMPLE_TOOLS]


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
    active_comp = [t for t in tools if t in PER_COMPARISON_TOOLS]
    active_persample = [t for t in tools if t in PER_SAMPLE_TOOLS]
    inactive = [t for t in config["tools"] if not config["tools"][t]["activate"]]

    logger.info("=" * 60)
    logger.info("NanoRNAMod — Active tools summary")
    logger.info(f"  Comparison tools ({len(active_comp)}): {', '.join(active_comp) or 'none'}")
    logger.info(f"  Per-sample tools ({len(active_persample)}): {', '.join(active_persample) or 'none'}")
    logger.info(f"  Inactive tools   ({len(inactive)}): {', '.join(inactive) or 'none'}")
    logger.info(f"  Samples: {', '.join(samples.index)}")
    logger.info(f"  Comparisons: {', '.join(comparisons)}")
    if config.get("benchmark", {}).get("truth_set", ""):
        logger.info(f"  Benchmarking: enabled (truth_set={config['benchmark']['truth_set']})")
        logger.info(f"    Comparison tools to benchmark: {', '.join(active_comp) or 'none'}")
        logger.info(f"    Per-sample tools to benchmark: {', '.join(active_persample) or 'none'}")
    else:
        logger.info("  Benchmarking: disabled (no truth_set)")
    logger.info("=" * 60)

    final_output = [f"{RESULT_ROOT}/workflow_version.json"]
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

    has_gtf = os.path.exists(config["reference"]["transcriptome_gtf"])

    # Comparison-based tool outputs
    for tool in tools:
        if tool in PER_COMPARISON_TOOLS:
            final_output += expand(
                f"{RESULT_ROOT}/modifications/{tool}/{{comp}}/{tool}_results.tsv",
                comp=comparisons,
            )
            if has_gtf:
                final_output += expand(
                    f"{RESULT_ROOT}/modifications/{tool}/{{comp}}/{tool}_annotated_results.tsv",
                    comp=comparisons,
                )
        elif tool in PER_SAMPLE_TOOLS:
            final_output += expand(
                f"{RESULT_ROOT}/modifications/{tool}/{{sample}}/{tool}_results.tsv",
                sample=list(samples.index),
            )
            if has_gtf:
                final_output += expand(
                    f"{RESULT_ROOT}/modifications/{tool}/{{sample}}/{tool}_annotated_results.tsv",
                    sample=list(samples.index),
                )

    # Benchmarking outputs - resource metrics
    final_output += [f"{RESULT_ROOT}/benchmarks/resource_summary.tsv"]
    final_output += [f"{RESULT_ROOT}/benchmarks/resource_by_tool.tsv"]
    final_output += [f"{RESULT_ROOT}/benchmarks/resource_by_tool.pdf"]

    # Accuracy benchmarking outputs (requires truth_set)
    if config.get("benchmark", {}).get("truth_set", ""):
        active_comp_tools = get_active_comparison_tools()

        # Cross-tool summaries (coverage and per-tool metrics are built
        # automatically as dependencies of the aggregate rule)
        for agg_file in [
            "accuracy_summary.tsv", "accuracy_summary_overall.tsv",
            "accuracy_summary_by_comparison.tsv", "by_tool.tsv",
            "best_scores.tsv",
        ]:
            final_output += [f"{RESULT_ROOT}/benchmarks/cross_tool/{agg_file}"]

        # PDF report (intermediate files like statistics/, sensitivity/,
        # coverage/, threshold_* are built automatically as dependencies
        # of R figures and the PDF report)
        final_output += [f"{RESULT_ROOT}/benchmarks/benchmark_report.pdf"]

        # Per-tool comparison figures (one representative file per tool×comp)
        final_output += expand(
            f"{RESULT_ROOT}/benchmarks/per_tool/{{tool}}/{{comp}}/figures/roc_curve.pdf",
            tool=active_comp_tools, comp=comparisons,
        )
        # Per-tool summary figures (one representative file per tool)
        final_output += expand(
            f"{RESULT_ROOT}/benchmarks/per_tool/{{tool}}/figures/accuracy_by_comparison.pdf",
            tool=active_comp_tools,
        )
        # Per-sample figures
        active_ps_tools = get_active_per_sample_tools()
        if active_ps_tools:
            final_output += expand(
                f"{RESULT_ROOT}/benchmarks/per_tool/{{tool}}/{{sample}}/figures/score_distribution.pdf",
                tool=active_ps_tools, sample=list(samples.index),
            )

        # Cross-tool comparison curves (overlaid ROC/PR, bar chart, score dist)
        for ct_fig in [
            "cross_tool_roc", "cross_tool_pr",
            "cross_tool_auroc_auprc", "cross_tool_score_distribution",
        ]:
            final_output += [f"{RESULT_ROOT}/benchmarks/cross_tool/figures/{ct_fig}.pdf"]
            final_output += [f"{RESULT_ROOT}/benchmarks/cross_tool/figures/{ct_fig}_data.tsv"]

        # Publication figures (flat layout with co-located data)
        main_fig_names = [
            "overall_accuracy", "precision_recall", "auroc",
            "native_vs_fair", "best_score", "coverage_sensitivity",
            "resource_usage", "tool_ranking",
        ]
        for i, name in enumerate(main_fig_names, 1):
            final_output += [f"{RESULT_ROOT}/benchmarks/figures/fig{i}_{name}.pdf"]
            final_output += [f"{RESULT_ROOT}/benchmarks/figures/fig{i}_{name}_data.tsv"]
        supp_fig_names = [
            "per_comparison", "native_vs_fair_detail",
            "threshold_robustness", "effect_sizes", "score_heatmap",
        ]
        for i, name in enumerate(supp_fig_names, 1):
            final_output += [f"{RESULT_ROOT}/benchmarks/figures/sfig{i}_{name}.pdf"]

    logger.debug("=" * 60)
    logger.debug(f"NanoRNAMod — Final output files ({len(final_output)} total)")
    logger.debug("=" * 60)
    for f in final_output:
        logger.debug(f"  {f}")
    logger.debug("=" * 60)

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

