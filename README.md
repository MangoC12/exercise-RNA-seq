# RNA sequencing exercise

## Prerequisites

- [fastp](https://github.com/OpenGene/fastp)
- [STAR](https://github.com/alexdobin/STAR)
- [RSEM](https://github.com/deweylab/RSEM)

## Preparation

1. Install prerequisites.
2. Download the reference sequences in `ref`.

## Instructions

1. Write a download script to obtain the paired-end RNA-seq data for
   the following samples from
   Sequence Read Archive [PRJNA681149](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA681149):
    - SUM149PT Parent clone, DMSO treated
    - SUM149PT C2 clone, DMSO treated

2. Build a reproducible pipeline to trim the reads with `fastp`,
   align the trimmmed reads with `STAR`, and
   quantify the expression profiles with `rsem`.

3. Describe key steps that need to be added to the pipeline so that 
   it can used in production.

