#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=logs/sra_%A_%a.out
#SBATCH --error=logs/sra_%A_%a.err
#SBATCH --array=1-126%14
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

ACCESSIONS="sra_accessions_final.txt"
SRA_DIR="sra_downloads"
FASTQ_DIR="fastq_folder"

mkdir -p "$SRA_DIR" "$FASTQ_DIR" logs

# Get this task's line
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSIONS")
DGRP=$(echo "$LINE" | cut -d',' -f1)
ACC=$(echo "$LINE"  | cut -d',' -f2 | tr -d '[:space:]')

echo "[INFO] Task $SLURM_ARRAY_TASK_ID: $DGRP - $ACC"

# Skip if already successfully downloaded
if ls "$FASTQ_DIR"/${DGRP}_*.fastq 1>/dev/null 2>&1; then
    echo "[SKIP] $DGRP - $ACC already exists"
    exit 0
fi

# Download into its own subdirectory to avoid race conditions
prefetch "$ACC" --output-directory "$SRA_DIR/$ACC"

# Find the SRA file ONLY within this accession's own directory
SRR_DIR=$(find "$SRA_DIR/$ACC" -name "*.sra" -o -name "*.sralite" | grep -v vdbcache | head -1 | xargs dirname 2>/dev/null)

if [ -z "$SRR_DIR" ]; then
    echo "[FAILED] prefetch did not produce an .sra file for $DGRP - $ACC"
    exit 1
fi

SRR=$(basename "$SRR_DIR")
echo "[INFO] Resolved $ACC -> $SRR"

# Extract fastq into a temp per-accession folder to avoid collisions
TMPFASTQ="$FASTQ_DIR/tmp_${ACC}"
mkdir -p "$TMPFASTQ"

fasterq-dump "$SRR_DIR" \
    --outdir "$TMPFASTQ" \
    --threads 4 \
    --progress

# Check fastq files were actually created
FASTQS=("$TMPFASTQ"/*.fastq)
if [ ${#FASTQS[@]} -eq 0 ] || [ ! -f "${FASTQS[0]}" ]; then
    echo "[FAILED] fasterq-dump produced no fastq files for $DGRP - $ACC"
    rm -rf "$TMPFASTQ"
    exit 1
fi

# Verify this SRR is not identical to an already-existing one
for f in "${FASTQS[@]}"; do
    FSIZE=$(stat -c%s "$f")
    FNAME=$(basename "$f")
    # Check if any existing fastq has the same size AND different DGRP name (potential duplicate)
    DUPE=$(find "$FASTQ_DIR" -maxdepth 1 -name "DGRP_*_${SRR}*.fastq" ! -path "$TMPFASTQ/*" 2>/dev/null | head -1)
    if [ -n "$DUPE" ]; then
        echo "[WARNING] $DGRP - $ACC ($SRR) appears to duplicate existing file: $DUPE"
        echo "[WARNING] This may indicate multiple SRX accessions map to the same SRR"
    fi
done

# Rename and move to final fastq folder: DGRP_21_SRR933566_1.fastq
for f in "${FASTQS[@]}"; do
    base=$(basename "$f")
    mv "$f" "$FASTQ_DIR/${DGRP}_${base}"
done

rm -rf "$TMPFASTQ"

# Final check
if ! ls "$FASTQ_DIR"/${DGRP}_*.fastq 1>/dev/null 2>&1; then
    echo "[FAILED] No fastq found after rename for $DGRP - $ACC"
    exit 1
fi

echo "[DONE] $DGRP - $ACC -> $SRR"

