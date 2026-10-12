#!/bin/bash
DIR="/mnt/d/xQTL_2025_Data/DGRP_XQTL/trimmed_fastqs"

for R1 in "$DIR"/*_merged_R1_001.fastq.gz; do
    BASE=$(basename "$R1" _merged_R1_001.fastq.gz)
    R1_FILE="$DIR/${BASE}_merged_R1_001.fastq.gz"
    R2_FILE="$DIR/${BASE}_merged_R2_001.fastq.gz"

    if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
        echo "Trimming: $BASE"
        trim_galore --paired --length 40 --max_n 1 -q 20 -j 8 \
            -o "/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/RawProcessing" \
            "$R1_FILE" "$R2_FILE"
    else
        echo "Missing file(s) for $BASE"
    fi
done

