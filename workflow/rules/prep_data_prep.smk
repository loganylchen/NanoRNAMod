localrules:
    link_fastq,
    link_blow5,


rule link_fastq:
    input:
        fastq=get_raw_fastq,
    output:
        fastq=temp("{project}/results/fastq/{sample}.fq.gz"),
    log:
        "logs/{project}/link_fastq/{sample}.log",
    threads: 1
    resources:
        mem_mb=1024 * 2,
    priority: 100
    shell:
        "exec >{log} 2>&1; "
        'src_md5=$(md5sum "{input.fastq}" | awk \'{{print $1}}\'); '
        'echo "src md5: $src_md5  {input.fastq}"; '
        'cp -f "{input.fastq}" "{output.fastq}"; '
        'dst_md5=$(md5sum "{output.fastq}" | awk \'{{print $1}}\'); '
        'echo "dst md5: $dst_md5  {output.fastq}"; '
        'if [ "$src_md5" != "$dst_md5" ]; then '
        '  echo "ERROR: md5 mismatch after copy" >&2; '
        '  rm -f "{output.fastq}"; '
        '  exit 1; '
        'fi; '
        "echo `date` "


rule link_blow5:
    input:
        blow5=get_raw_blow5,
    output:
        blow5=temp("{project}/results/blow5/{sample}.blow5"),
    log:
        "logs/{project}/link_blow5/{sample}.log",
    threads: 1
    resources:
        mem_mb=1024 * 2,
    priority: 100
    shell:
        "exec >{log} 2>&1; "
        'src_md5=$(md5sum "{input.blow5}" | awk \'{{print $1}}\'); '
        'echo "src md5: $src_md5  {input.blow5}"; '
        'cp -f "{input.blow5}" "{output.blow5}"; '
        'dst_md5=$(md5sum "{output.blow5}" | awk \'{{print $1}}\'); '
        'echo "dst md5: $dst_md5  {output.blow5}"; '
        'if [ "$src_md5" != "$dst_md5" ]; then '
        '  echo "ERROR: md5 mismatch after copy" >&2; '
        '  rm -f "{output.blow5}"; '
        '  exit 1; '
        'fi; '
        "echo `date` "
