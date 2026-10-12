#!/bin/bash
INPUT_DIR="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/RawProcessing/bams"
OUTPUT_DIR="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/RawProcessing/bams/DuplicateClear"

mkdir -p "$OUTPUT_DIR"

for BAM in "$INPUT_DIR"/*.bam; do
    BASENAME=$(basename "$BAM" .bam)
    FINAL_OUTPUT="$OUTPUT_DIR/${BASENAME}_dedup.bam"

    if [[ -f "$FINAL_OUTPUT" ]]; then
        echo "✅ Already processed: $FINAL_OUTPUT"
        continue
    fi

    echo "🔄 Processing $BAM..."

    samtools sort -n -o "$OUTPUT_DIR/${BASENAME}_namesorted.bam" "$BAM"
    samtools fixmate -m "$OUTPUT_DIR/${BASENAME}_namesorted.bam" "$OUTPUT_DIR/${BASENAME}_fixmate.bam"
    samtools sort -o "$OUTPUT_DIR/${BASENAME}_positionsorted.bam" "$OUTPUT_DIR/${BASENAME}_fixmate.bam"
    samtools markdup -r -s "$OUTPUT_DIR/${BASENAME}_positionsorted.bam" "$FINAL_OUTPUT" 2> "$OUTPUT_DIR/${BASENAME}_dedup_stats.txt"

    echo "✅ Output written to $FINAL_OUTPUT"

    rm "$OUTPUT_DIR/${BASENAME}_namesorted.bam"
    rm "$OUTPUT_DIR/${BASENAME}_fixmate.bam"
    rm "$OUTPUT_DIR/${BASENAME}_positionsorted.bam"
done
