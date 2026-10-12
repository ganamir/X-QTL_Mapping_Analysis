#!/bin/bash
bam_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/bams"
out_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams"
mkdir -p "$out_dir"

# Get unique line numbers
ls "$bam_dir"/*_aligned.bam | grep -oP 'DGRP_\d+' | sort -u | while read line; do
    bams=("$bam_dir"/${line}_*_aligned.bam)
    deduped_bams=()

    # Dedup each SRR BAM separately
    for bam in "${bams[@]}"; do
        base=$(basename "$bam" .bam)
        deduped="$out_dir/${base}_dedup.bam"
        deduped_bams+=("$deduped")

        samtools sort -n -o "$out_dir/${base}_namesorted.bam" "$bam"
        samtools fixmate -m "$out_dir/${base}_namesorted.bam" "$out_dir/${base}_fixmate.bam"
        samtools sort -o "$out_dir/${base}_positionsorted.bam" "$out_dir/${base}_fixmate.bam"
        samtools markdup -r -s "$out_dir/${base}_positionsorted.bam" "$deduped" 2> "$out_dir/${base}_dedup_stats.txt"

        rm "$out_dir/${base}_namesorted.bam" "$out_dir/${base}_fixmate.bam" "$out_dir/${base}_positionsorted.bam"
    done

    # Merge deduped BAMs per line
    final="$out_dir/${line}.bam"
    if [ "${#deduped_bams[@]}" -gt 1 ]; then
        echo "Merging ${#deduped_bams[@]} deduped BAMs for $line"
        samtools merge -f "$final" "${deduped_bams[@]}"
        rm "${deduped_bams[@]}"
    else
        mv "${deduped_bams[0]}" "$final"
    fi

    samtools index "$final"
    echo "Done: $final"
done

