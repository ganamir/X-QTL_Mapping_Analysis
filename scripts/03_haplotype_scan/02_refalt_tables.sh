#!/bin/bash
#SBATCH --job-name=refalt
#SBATCH --output=/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups/haplotype_inf_final/process/logs/refalt_%a.out
#SBATCH --array=1-5
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G

ref="/mnt/d/xQTL_2025_Data/ref/dm6.fa"
bam_dir_founders="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams_RG"
bam_dir_pool="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups"
output="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups/haplotype_inf_final/process"
mkdir -p "$output/logs"

ls "$bam_dir_founders"/DGRP_*.bam "$bam_dir_pool"/{A,C}*.bam > /tmp/all_bams.txt

chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "Processing $mychr"

bcf_out="${output}/calls.${mychr}.bcf"

bcftools mpileup -I -d 1000 -r "$mychr" -a "FORMAT/AD,FORMAT/DP" -f "$ref" -b /tmp/all_bams.txt \
    --threads 10 | \
    bcftools call -mv -Ob -o "$bcf_out"

bcftools index "$bcf_out"

echo -ne "CHROM\tPOS" > "${output}/RefAlt.${mychr}.txt"
bcftools query -l "$bcf_out" | awk '{printf("\tREF_%s\tALT_%s",$1,$1)}' >> "${output}/RefAlt.${mychr}.txt"
echo -ne "\n" >> "${output}/RefAlt.${mychr}.txt"

bcftools view -m2 -M2 -v snps -i 'QUAL>59' "$bcf_out" | \
    bcftools query -f '%CHROM\t%POS[\t%AD{0}\t%AD{1}]\n' | \
    grep -v '\.' >> "${output}/RefAlt.${mychr}.txt"

echo "Finished $mychr"

