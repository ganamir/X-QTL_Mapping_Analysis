#!/bin/bash
#SBATCH --job-name=add_RG
#SBATCH --output=/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/logs/RG_%a.out
#SBATCH --array=1-100%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

in_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams"
out_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams_RG"
mkdir -p "$out_dir"

# Get the Nth BAM for this array task
bam=$(ls "$in_dir"/DGRP_*.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
line=$(basename "$bam" .bam)
line_num=$(echo "$line" | cut -d_ -f2)

picard AddOrReplaceReadGroups \
    I="$bam" \
    O="$out_dir/${line}.bam" \
    RGID=${line_num} \
    RGSM=DGRP_${line_num} \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGLB=Lib1 \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT

samtools index "$out_dir/${line}.bam"
echo "Done: $line"

