rule epinano:
    input:
        control_bam="results/alignments/{control}_filtered.bam",
        control_bai="results/alignments/{control}_filtered.bam.bai",
        native_bam="results/alignments/{native}_filtered.bam",
        native_bai="results/alignments/{native}_filtered.bam.bai",
        reference=config['reference']['transcriptome_fasta'],
        reference_dict=config['reference']['transcriptome_fasta']+'.dict'
    output:
        "results/eligos2/{native}_{control}"
    params:
        prefix="{native}_{control}",
        extra=config['params']['eligo2']
    threads: config['threads']['eligo2']
    log:
        stdout="logs/eligos2/{native}_{control}.log"

    benchmark:
        "benchmarks/{native}_{control}.eligos2.benchmark.txt"
    container:
        "docker://btrspg/eligos2:latest"
    shell:
        "eligos2 pair_diff_mod -tbam {input.native_bam} "
        "-cbam {input.control_bam} "
        "-ref {input.reference} "
        "-t {threads} "
        "-o {output} {params.extra}"

rule epinano_prep:
    input:
        sample_bam="results/alignments/{sample}_filtered.bam",
        sample_bai="results/alignments/{sample}_filtered.bam.bai",
        reference=config['reference']['transcriptome_fasta'],
        reference_dict=config['reference']['transcriptome_fasta'] + '.dict'
    output:
        "results/dataprep/{sample}_epinano_dataprep/"
    params:
        extra=config['params']['epinano_dataprep']
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
        "-T t "


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