#!/usr/bin/env bash
# ============================================================
# run_eggnog.sh — 对蛋白 FASTA 运行 eggNOG-mapper（离线注释）
#
# 依赖：
#   * conda 环境（默认名 eggnog）中安装 eggnog-mapper 与 diamond
#   * eggNOG 数据库目录，含 eggnog_proteins.dmnd
#
# 环境变量（均可覆盖）：
#   PROTEIN_FASTA  蛋白 FASTA
#   EGGNOG_DATA    eggNOG 数据库目录
#   OUTDIR         输出目录
#   EGGNOG_ENV     conda 环境名
#   EGGNOG_CPU     线程数
#
# 例：
#   PROTEIN_FASTA=/path/to/protein.faa \
#   EGGNOG_DATA=~/eggnog_data \
#   EGGNOG_CPU=8 \
#   bash Pipeline/run_eggnog.sh
# ============================================================
set -euo pipefail

PROTEIN_FASTA="${PROTEIN_FASTA:-$HOME/rna-seq/data/reference/protein.faa}"
EGGNOG_DATA="${EGGNOG_DATA:-$HOME/eggnog_data}"
OUTDIR="${OUTDIR:-$EGGNOG_DATA/emapper}"
EGGNOG_ENV="${EGGNOG_ENV:-eggnog}"
EGGNOG_CPU="${EGGNOG_CPU:-8}"

[ -f "$PROTEIN_FASTA" ] || { echo "找不到蛋白 FASTA: $PROTEIN_FASTA" >&2; exit 1; }
[ -f "$EGGNOG_DATA/eggnog_proteins.dmnd" ] || {
  echo "缺少 eggNOG 数据库：$EGGNOG_DATA/eggnog_proteins.dmnd" >&2
  echo "请先运行 download_eggnog_data.py" >&2
  exit 1
}

mkdir -p "$OUTDIR"
echo "[eggnog] emapper.py -i $PROTEIN_FASTA --cpu $EGGNOG_CPU"
conda run --no-capture-output -n "$EGGNOG_ENV" emapper.py \
  -i "$PROTEIN_FASTA" \
  -o "$OUTDIR/out" \
  --cpu "$EGGNOG_CPU" \
  --data_dir "$EGGNOG_DATA"

echo "[done] $OUTDIR/out.emapper.annotations"
