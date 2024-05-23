rule epinano_prep:
    input:
        sample_bam="results/alignments/{sample}_filtered.bam",
        sample_bai="results/alignments/{sample}_filtered.bam.bai",
        reference=config['reference']['transcriptome_fasta'],
        reference_dict=config['reference']['transcriptome_fasta'] + '.dict'
    output:
        per_site = "results/alignments/{sample}_filtered.plus_strand.per.site.csv",
        kmer_5_site = "results/alignments/{sample}_filtered.plus_strand.per.site.5mer.csv",
    params:
        extra=config['params']['epinano_dataprep'],
        prefix="results/alignments/{sample}_filtered"
    threads: config['threads']['epinano']
    log:
        "logs/epinano_prep/{sample}.log"
    benchmark:
        "benchmarks/{sample}.epinano_prep.benchmark.txt"
    container:
        "docker://btrspg/epinano:latest"
    shell:
        "epinano_variants -R {input.reference} "
        "-b {input.sample_bam} "
        "-n {threads} "
        "-T t 2>{log} && "
        "python3 /opt/EpiNano/misc/Slide_Variants.py {output.per_site} 5 2>>{log} &&"
        "python3 /opt/EpiNano/Epinano_Predict.py "
        "--model /opt/EpiNano/models/rrach.q3.mis3.del3.linear.dump "
        "--predict {output.kmer_5_site} "
        "--columns 8,13,23 "
        "--out_prefix {params.prefix}  2>>{log}"

# rule epinano_

rule epinano_prep_sampled:
    input:
        sample_bam="results/alignments/{sample}_filtered_{sample_size}_{n}.bam",
        sample_bai="results/alignments/{sample}_filtered_{sample_size}_{n}.bam.bai",
        reference=config['reference']['transcriptome_fasta'],
        reference_dict=config['reference']['transcriptome_fasta'] + '.dict'
    output:
        per_site = "results/alignments/{sample}_filtered_{sample_size}_{n}.plus_strand.per.site.csv",
        kmer_5_site = "results/alignments/{sample}_filtered_{sample_size}_{n}.plus_strand.per.site.5mer.csv",
    params:
        extra=config['params']['epinano_dataprep'],
        prefix="results/alignments/{sample}_filtered_{sample_size}_{n}"
    threads: config['threads']['epinano']
    log:
        "logs/epinano_prep/{sample}_{sample_size}_{n}.log"
    benchmark:
        "benchmarks/{sample}_{sample_size}_{n}.epinano_prep.benchmark.txt"
    container:
        "docker://btrspg/epinano:latest"
    shell:
        "epinano_variants -R {input.reference} "
        "-b {input.sample_bam} "
        "-n {threads} "
        "-T t 2>{log} && "
        "python3 /opt/EpiNano/misc/Slide_Variants.py {output.per_site} 5 2>>{log} &&"
        "python3 /opt/EpiNano/Epinano_Predict.py "
        "--model /opt/EpiNano/models/rrach.q3.mis3.del3.linear.dump "
        "--predict {output.kmer_5_site} "
        "--columns 8,13,23 "
        "--out_prefix {params.prefix}  2>>{log}"




# rule epinano:
#     input:
#         control="results/alignments/{control}_filtered.plus_strand.per.site.csv",
#         native="results/alignments/{native}_filtered.plus_strand.per.site.csv",
#     output:
#         "results/alignments/{wt}_{ko}_delta.5mer.csv"
#     params:
#         extra=config['params']['epinano']
#     threads: config['threads']['epinano']
#     log:
#         "logs/epinano/{wt}_{ko}.log"
#     benchmark:
#         "benchmarks/{wt}_{ko}.epinano.benchmark.txt"
#     container:
#         "docker://btrspg/epinano:latest"
#     shell:
#         "epinano_predict --model {input.model} "
#         "--predict {input.delta} "
#         "--columns 7,12,22 "
#         "--out_prefix {output} 2>{log}"

'''Epinano=path_to_Epinano_software
ref=path_to_reference_transcriptome
wt_bam=path_to_bam_file
ko_bam=path_to_bam_file

#extract features
python $Epinano/Epinano_Variants.py -R $ref \
-b $wt_bam \
-s $Epinano/misc/sam2tsv.jar \
-n 16 \
-T t

#slide features
python $Epinano/misc/Slide_Variants.py wt.plus_strand.per.site.csv 5
# output results/alignments/RNA081120181_filtered.plus_strand.per.site.5mer.csv



#predict modifications
python $Epinano/Epinano_Predict.py \
--model $Epinano/models/rrach.q3.mis3.del3.linear.dump \
--predict wt.plus_strand.per.site.5mer.csv \
--columns 8,13,23 \
--out_prefix wt

#apply the same analysis in ko


python $Epinano/misc/Epinano_make_delta.py $epinano/wt.plus_strand.per.site.5mer.csv \
$epinano/ko.plus_strand.per.site.5mer.csv \
5 5 > wt_ko_delta.5mer.csv

#predict modifications
python $Epinano/Epinano_Predict.py \
--model $Epinano/models/rrach.deltaQ3.deltaMis3.deltaDel3.linear.dump \
--predict wt_ko_delta.5mer.csv \
--columns 7,12,22 \
--out_prefix com'''