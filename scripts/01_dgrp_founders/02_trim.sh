#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --output=logs/trim_%A_%a.out
#SBATCH --error=logs/trim_%A_%a.err
#SBATCH --array=1-100%5
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

FASTQ_DIR="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/fastq_folder"
TRIM_DIR="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/trimmed_fastqs"

mkdir -p "$TRIM_DIR" logs

DGRP=$(ls "$FASTQ_DIR"/DGRP_*.fastq | grep -oP 'DGRP_\d+' | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "[INFO] Task $SLURM_ARRAY_TASK_ID: $DGRP"

# Skip if already trimmed
if ls "$TRIM_DIR"/${DGRP}_*trimmed* 1>/dev/null 2>&1; then
    echo "[SKIP] $DGRP already trimmed"
    exit 0
fi

# Paired-end runs
for R1 in "$FASTQ_DIR"/${DGRP}_*_1.fastq; do
    [ -f "$R1" ] || continue
    R2="${R1/_1.fastq/_2.fastq}"
    if [[ -f "$R2" ]]; then
        echo "[PE] Trimming $R1 + $R2"
        trim_galore --paired --length 30 --max_n 1 -q 20 -j 8 \
            -o "$TRIM_DIR" "$R1" "$R2"
    fi
done

# Single-end runs
for SE in "$FASTQ_DIR"/${DGRP}_*.fastq; do
    [[ "$SE" == *_1.fastq ]] && continue
    [[ "$SE" == *_2.fastq ]] && continue
    echo "[SE] Trimming $SE"
    trim_galore --length 30 --max_n 1 -q 20 -j 8 \
        -o "$TRIM_DIR" "$SE"
done

echo "[DONE] $DGRP"


