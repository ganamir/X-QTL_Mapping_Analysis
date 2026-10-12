#!/bin/bash

# Run directly with: bash merge_lanes_simple.sh

SOURCE_DIR="/mnt/d/xQTL_2025_Data/DGRP_XQTL"
OUTPUT_DIR="/mnt/d/xQTL_2025_Data/DGRP_XQTL/trimmed_fastqs"

echo "=== Starting Lane Merge ==="
echo "Source: $SOURCE_DIR"
echo "Output: $OUTPUT_DIR"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Counter
TOTAL=0
SUCCESS=0
FAILED=0

# Loop through all unique samples
for L001_R1 in "$SOURCE_DIR"/2_[A-E]*_L001_R1_001.fastq.gz; do
    
    # Extract base name
    BASE=$(basename "$L001_R1")
    BASE=${BASE%%_L001_R1_001.fastq.gz}
    
    ((TOTAL++))
    
    echo "[$TOTAL] Processing: $BASE"
    
    # Define file paths
    L001_R1="$SOURCE_DIR/${BASE}_L001_R1_001.fastq.gz"
    L001_R2="$SOURCE_DIR/${BASE}_L001_R2_001.fastq.gz"
    L002_R1="$SOURCE_DIR/${BASE}_L002_R1_001.fastq.gz"
    L002_R2="$SOURCE_DIR/${BASE}_L002_R2_001.fastq.gz"
    
    # Check all files exist
    if [[ -f "$L001_R1" && -f "$L001_R2" && -f "$L002_R1" && -f "$L002_R2" ]]; then
        
        # Merge R1
        cat "$L001_R1" "$L002_R1" > "$OUTPUT_DIR/${BASE}_merged_R1_001.fastq.gz"
        
        # Merge R2
        cat "$L001_R2" "$L002_R2" > "$OUTPUT_DIR/${BASE}_merged_R2_001.fastq.gz"
        
        echo "  ✓ Complete"
        ((SUCCESS++))
    else
        echo "  ✗ ERROR: Missing files"
        ((FAILED++))
    fi
done

echo ""
echo "=== Summary ==="
echo "Total samples: $TOTAL"
echo "Successful: $SUCCESS"
echo "Failed: $FAILED"
echo ""
echo "Output location: $OUTPUT_DIR"
ls -lh "$OUTPUT_DIR"/*_merged_*.fastq.gz 2>/dev/null | wc -l | xargs echo "Merged files created:"

