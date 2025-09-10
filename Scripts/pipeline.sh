#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#activate conda environment
#conda init
#conda activate RNAseq_env

# 解压出基因组 FASTA
if [ ! -s ref/genome.fa ]; then
    gunzip -c ref/GRCh38.primary_assembly.genome.fa.gz > ref/genome.fa
else
    echo "skipping"
fi

if [ ! -s ref/transcripts.fa ]; then
    gunzip -c ref/gencode.v42.transcripts.fa.gz > ref/transcripts.fa
else
    echo "skipping"
fi

if [ ! -s ref/annotation.gtf ]; then
    gunzip -c ref/gencode.v42.annotation.gtf.gz > ref/annotation.gtf
else
    echo "skipping"
fi

# Reference files
GENOME_FA="ref/genome.fa"
GTF="ref/annotation.gtf" 
GENOME_NAME="genome"  # 注:RSEM 索引名不要带路径
THREADS=6

# Define directories
RAW_DIR="fq"
TRIMMED_DIR="ref/trimmed" #fastp输出/指控后数据
STAR_INDEX="ref/star_index" #STAR索引
RSEM_INDEX="ref/rsem_index" #RSEM索引
QUANT_DIR="ref/quant" #RSEM定量输出

mkdir -p "$TRIMMED_DIR" "$STAR_INDEX" "$RSEM_INDEX" "$QUANT_DIR"



# 1. 构建 STAR 索引
if [ ! -s ${STAR_INDEX}/Genome ]; then
  echo "[STAR] Building genome index..."
  STAR --runThreadN ${THREADS} \
       --runMode genomeGenerate \
       --genomeDir "$STAR_INDEX" \
       --genomeFastaFiles "$GENOME_FA" \
       --sjdbGTFfile "$GTF" \
       --sjdbOverhang 100 \
       --limitGenomeGenerateRAM 24000000000 \
       --genomeSAindexNbases 12
else
  echo "[STAR] Index exists. Skipping."
fi

# 2. 构建 RSEM 索引
if [ ! -s ${RSEM_INDEX}/${GENOME_NAME}.grp ]; then
  echo "[RSEM] Building reference..."
  rsem-prepare-reference \
    --gtf "$GTF" \
    --star \
    --star-path "$(dirname $(which STAR))" \
    "$GENOME_FA" "$RSEM_INDEX/$GENOME_NAME"
else
  echo "[RSEM] Reference exists. Skipping."
fi

# 3. 构建样本SAMPLES列表
# 兼容 *_1.fastq 与 *_1.fastq.gz，两者都支持；输出统一用 .fastq.gz
mapfile -t SAMPLES < <(ls ${RAW_DIR}/*_1.fastq* 2>/dev/null | sed 's#.*/##; s/_1\.fastq.*$//')
if [ ${#SAMPLES[@]} -eq 0 ]; then
  echo "ERROR: can't find ${RAW_DIR}/ *_1.fastq 或 *_1.fastq.gz" >&2
  exit 1
fi

# 4. 循环：fastp 修剪 → RSEM+STAR 对齐定量
for SAMPLE in \"${SAMPLES[@]}\"; 
do
  echo "=== Processing sample: $SAMPLE ==="

  # 输入 R1/R2（可能是 .fastq 或 .fastq.gz）
  R1_GZ="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
  R2_GZ="${RAW_DIR}/${SAMPLE}_2.fastq.gz"
  R1_FQ="${RAW_DIR}/${SAMPLE}_1.fastq"
  R2_FQ="${RAW_DIR}/${SAMPLE}_2.fastq"

  if   [ -s "${R1_GZ}" ] && [ -s "${R2_GZ}" ]; then
    IN1="${R1_GZ}"; IN2="${R2_GZ}"
  elif [ -s "${R1_FQ}" ] && [ -s "${R2_FQ}" ]; then
    IN1="${R1_FQ}"; IN2="${R2_FQ}"
  else
    echo "WARN: Can't find ${SAMPLE} 's two sides reads, skipping." >&2
    continue
  fi

  # fastp 输出统一为 .fastq.gz
  OUT1="${TRIMMED_DIR}/${SAMPLE}_1.trimmed.fastq.gz"
  OUT2="${TRIMMED_DIR}/${SAMPLE}_2.trimmed.fastq.gz"

  # 4.1 fastp（若已存在结果则跳过）
  if [ ! -s "${OUT1}" ] || [ ! -s "${OUT2}" ]; then
    echo "[fastp] ${SAMPLE}"
    fastp \
      -i "${IN1}" -I "${IN2}" \
      -o "${OUT1}" -O "${OUT2}" \
      --thread "${THREADS}" \
      --compression 6 \
      --html "${TRIMMED_DIR}/${SAMPLE}.fastp.html" \
      --json "${TRIMMED_DIR}/${SAMPLE}.fastp.json"
  else
    echo "[fastp] ${SAMPLE} trimmed files exist. Skipping."
  fi

  # 4.2 RSEM+STAR：RSEM 内部调用 STAR 完成对齐与定量（若已完成则跳过）
  # 说明：--star-gzipped-read-file 告知 RSEM 输入是 .gz
  if [ ! -s "${QUANT_DIR}/${SAMPLE}.genes.results" ]; then
    echo "[RSEM+STAR] ${SAMPLE}"
    rsem-calculate-expression \
      --paired-end \
      --star \
      --star-path "$(dirname "$(which STAR)")" \
      --star-gzipped-read-file \
      --output-genome-bam \
      -p "${THREADS}" \
      "${OUT1}" "${OUT2}" \
      "${RSEM_INDEX}/genome" \
      "${QUANT_DIR}/${SAMPLE}"
  else
    echo "[RSEM] ${SAMPLE} results exist. Skipping."
  fi
done

# 5. 结果汇总（计数矩阵 & TPM矩阵）
# RSEM 自带汇总工具，expected_count 用 rsem-generate-data-matrix 直接汇总
echo "[RSEM] Generate count matrix..."
mapfile -t GENE_RESULTS < <(ls "${QUANT_DIR}"/*.genes.results 2>/dev/null | sort)
if [ ${#GENE_RESULTS[@]} -gt 0 ]; then
  rsem-generate-data-matrix "${GENE_RESULTS[@]}" > "${QUANT_DIR}/genes.expected_count.matrix.txt"
else
  echo "WARN: Can't find *.genes.results, skipping." >&2
fi

# TPM 汇总（rsem-generate-data-matrix 没有直接的 TPM 选项；这里用简单 awk 合并列6）
# 第一列为基因ID（来自第一个样本），后续每列是对应样本TPM
echo "[RSEM] Generate TPM matrix..."
if [ ${#GENE_RESULTS[@]} -gt 0 ]; then
  # 取样本名（不带路径与后缀）
  SAMPLES_PRINT=()
  for f in "${GENE_RESULTS[@]}"; do
    base=$(basename "$f")
    SAMPLES_PRINT+=("${base%.genes.results}")
  done

  # 生成表头
  {
    printf "gene_id"
    for s in "${SAMPLES_PRINT[@]}"; do printf "\t%s" "$s"; done
    printf "\n"
  } > "${QUANT_DIR}/genes.tpm.matrix.txt"

  # 用第一个文件的 gene_id 当行框架
  awk 'NR>1{print $1}' "${GENE_RESULTS[0]}" > "${QUANT_DIR}/._gene_ids.tmp"

  # 逐个样本抽 TPM（第6列）并拼到一起
  paste_files=("${QUANT_DIR}/._gene_ids.tmp")
  for f in "${GENE_RESULTS[@]}"; do
    awk 'NR>1{print $6}' "$f" > "${f}.tpm.col"
    paste_files+=("${f}.tpm.col")
  done

  paste "${paste_files[@]}" >> "${QUANT_DIR}/genes.tpm.matrix.txt"
  rm -f "${QUANT_DIR}/._gene_ids.tmp" "${QUANT_DIR}"/*.tpm.col
else
  echo "WARN: Can't find *.genes.results, skipping." >&2
fi

echo "RNA-seq pipeline finished. Results in: ${QUANT_DIR}"