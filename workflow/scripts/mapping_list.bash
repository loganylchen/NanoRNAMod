#! /bin/bash


helpFunction()
{
   echo ""
   echo "Usage: $0 -m mapping_list -d mapping_dir -t threads "
   echo -e "\t-m mapping list file"
   echo -e "\t-d mapping_dir"
   echo -e "\t-t threads"
   echo -e "\t-b bam file"
   echo -e "\t-l bam list file"
   echo -e "\t-p minimap2 parameters "
   exit 1 # Exit script after printing help
}

while getopts "m:d:t:h:b:l:p:r:" opt
do
   case "$opt" in
      m ) mapping_list="$OPTARG"
          echo "mapping_list: $mapping_list"
          ;;
      d ) mapping_dir="$OPTARG"
          echo "mapping_dir: $mapping_dir"
          ;;
      t ) threads="$OPTARG"
          echo "threads: $threads"
          ;;
      b ) output_bam="$OPTARG"
          echo "bam: $output_bam"
          ;;
      l ) bam_list="$OPTARG"
          echo "bam_list: $bam_list"
          ;;
      p ) parameter="$OPTARG"
          echo "parameter: $parameter"
          ;;
      r ) ref="$OPTARG"
          echo "ref: $ref"
          ;;
      h ) helpFunction ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


for line in $(cat $mapping_list)
do
    transcript=$(echo $line|awk -F "\t" "print $1");
    fastq=$(echo $line|awk -F "\t" "print $2");
    fasta=$(echo $line|awk -F "\t" "print $3");
    bam=$(echo $line|awk -F "\t" "print $4");
    echo "Working on $transcript" ;
    minimap2 -t $threads $parameter $fasta $fastq | samtools view -Sbh | samtools sort - -o $bam;
    samtools index $bam
done
awk -F "\t" "print $4" $mapping_list > $bam_list
awk -F "\t" "print @SQ\tSN:$1\tLN:$2" $ref.fai > $mapping_dir/header.tmp
samtools merge -b $bam_list -h $mapping_dir/header.tmp | samtools sort - -o $output_bam ;
samtools index $output_bam
rm $mapping_dir/header.tmp