FROM snakemake/snakemake:v7.31.0

RUN apt update && \
    apt install -y build-essential  libz-dev && \
    pip install pyfastx interlap mappy