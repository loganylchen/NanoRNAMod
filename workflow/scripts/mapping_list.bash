#!/usr/bin/env bash


exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

__help="
This is my first bash script for snakemake
So I want to print some unused sentences to test if the program is running well

Here are the inputs, outputs, parameters
=================================================================================

input:
    mapping_list=${snakemake_input[mapping_list]},
    mapping_dir=${snakemake_input[mapping_dir]},
output:
    bam_list=${snakemake_output[bam_list]},
    bam=${snakemake_output[bam]},
log: ${snakemake_log[0]}


=================================================================================

"

echo __help >> ${snakemake_log[0]}


for line in $(cat ${snakemake_input[mapping_list]})
do
    transcript=$(echo $line|awk -F "\t" "print $1");
    fastq=$(echo $line|awk -F "\t" "print $2");
    fasta=$(echo $line|awk -F "\t" "print $3");
    bam=$(echo $line|awk -F "\t" "print $4");
    echo "Working on $transcript" ;
    minimap2 -t ${snakemake[threads]} ${snakemake_params[extra]} $fasta $fastq | samtools view -Sbh | samtools sort - -o $bam;
    samtools index $bam
done
awk -F "\t" "print $4" ${snakemake_input[mapping_list]} > ${snakemake_output[bam_list]}
awk -F "\t" "print @SQ\tSN:$1\tLN:$2" ${snakemake_params[transcriptome_fasta]}.fai > ${snakemake_input[mapping_dir]}/header.tmp
samtools merge -b ${snakemake_output[bam_list]} -h ${snakemake_input[mapping_dir]}/header.tmp | samtools sort - -o ${snakemake_output[bam]} ;
samtools index ${snakemake_output[bam_list]}
rm ${snakemake_input[mapping_dir]}/header.tmp








