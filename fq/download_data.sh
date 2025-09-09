# TODO Download data from SRA

# (first time) conda install -c bioconda sra-tools -y # for prefetch, fasterq-dump

#!/bin/bash
# download SRA data & turn it into fastq
# Sequence Read Archive [PRJNA681149](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA681149):
    - SUM149PT Parent clone, DMSO treated #"SRR13155444" "SRR13155443" "SRR13155442"
    - SUM149PT C2 clone, DMSO treated #"SRR13155429" "SRR13155428" "SRR13155427"

conda activate RNAseq_env

SAMPLES=("SRR13155429" "SRR13155428" "SRR13155427" "SRR13155444" "SRR13155443" "SRR13155442")

for SAMPLE in "${SAMPLES[@]}"
do
    if [[ -f fq/${SAMPLE}_1.fastq && -f fq/${SAMPLE}_2.fastq ]]; then
        continue
    fi
    prefetch "$SAMPLE" > /dev/null || continue
    fasterq-dump --split-files --outdir fq "$SAMPLE" > /dev/null || continue
done