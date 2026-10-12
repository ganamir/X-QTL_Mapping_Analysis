#!/bin/bash
indir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3"
outRG="${indir}/RG_groups"
mkdir -p "$outRG"

# Treatment samples A1-A8 (A6-A12 + B1)
declare -A t_samples
t_samples=([A6]=S120 [A7]=S124 [A8]=S128 [A9]=S132 [A10]=S136 [A11]=S140 [A12]=S144 [B1]=S98)
rg_num=1
for sample in A6 A7 A8 A9 A10 A11 A12 B1; do
    s=${t_samples[$sample]}
    inbam="2_${sample}_${s}_merged_aligned_dedup.bam"
    shortname=$(basename "$inbam" .bam)
    rgid="A${rg_num}"
    rgsm="A${rg_num}"
    echo "Processing $inbam -> RGID=$rgid"
    picard AddOrReplaceReadGroups \
        I="${indir}/${inbam}" \
        O="${outRG}/${shortname}.RGfixed.bam" \
        SORT_ORDER=coordinate \
        RGPL=illumina RGPU=D109LACXX RGLB=Lib1 \
        RGID=${rgid} RGSM=${rgsm} \
        VALIDATION_STRINGENCY=LENIENT
    mv "${outRG}/${shortname}.RGfixed.bam" "${outRG}/A${rg_num}.bam"
    ((rg_num++))
done

# Control samples C1-C8 (B2-B9)
declare -A c_samples
c_samples=([B2]=S103 [B3]=S108 [B4]=S113 [B5]=S117 [B6]=S121 [B7]=S125 [B8]=S129 [B9]=S133)
rg_num=1
for sample in B2 B3 B4 B5 B6 B7 B8 B9; do
    s=${c_samples[$sample]}
    inbam="2_${sample}_${s}_merged_aligned_dedup.bam"
    shortname=$(basename "$inbam" .bam)
    rgid="C${rg_num}"
    rgsm="C${rg_num}"
    echo "Processing $inbam -> RGID=$rgid"
    picard AddOrReplaceReadGroups \
        I="${indir}/${inbam}" \
        O="${outRG}/${shortname}.RGfixed.bam" \
        SORT_ORDER=coordinate \
        RGPL=illumina RGPU=D109LACXX RGLB=Lib1 \
        RGID=${rgid} RGSM=${rgsm} \
        VALIDATION_STRINGENCY=LENIENT
    mv "${outRG}/${shortname}.RGfixed.bam" "${outRG}/C${rg_num}.bam"
    ((rg_num++))
done



# index the renamed BAMs
for b in "${outRG}"/*.bam; do samtools index "$b"; done
