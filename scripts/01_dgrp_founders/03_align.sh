input_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/trimmed_fastqs"
output_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/bams"
mkdir -p "$output_dir"
reference_genome="/mnt/d/xQTL_2025_Data/ref/dm6.fa"
threads=50

# --- Paired-end ---
for file1 in "$input_dir"/DGRP_*_1_val_1.fq; do
    [ -f "$file1" ] || continue
    file2="${file1/_1_val_1.fq/_2_val_2.fq}"
    if [ ! -f "$file2" ]; then
        echo "R2 not found for $file1, skipping."
        continue
    fi

    base=$(basename "$file1" _1_val_1.fq)
    line_num=$(echo "$base" | cut -d_ -f2)
    rg="@RG\tID:${line_num}\tSM:${base}\tPL:ILLUMINA"

    temp_bam="$output_dir/${base}_temp_aligned.bam"
    out_bam="$output_dir/${base}_aligned.bam"

    echo "PE: $file1 + $file2"
    bwa mem -M -t "$threads" -v 3 -R "$rg" "$reference_genome" "$file1" "$file2" \
        | samtools view -bS - > "$temp_bam"
    samtools sort "$temp_bam" -o "$out_bam"
    samtools index "$out_bam"
    rm "$temp_bam"
done

# --- Single-end ---
for file1 in "$input_dir"/DGRP_*_trimmed.fq; do
    [ -f "$file1" ] || continue

    base=$(basename "$file1" _trimmed.fq)
    line_num=$(echo "$base" | cut -d_ -f2)
    rg="@RG\tID:${line_num}\tSM:${base}\tPL:ILLUMINA"

    temp_bam="$output_dir/${base}_temp_aligned.bam"
    out_bam="$output_dir/${base}_aligned.bam"

    echo "SE: $file1"
    bwa mem -M -t "$threads" -v 3 -R "$rg" "$reference_genome" "$file1" \
        | samtools view -bS - > "$temp_bam"
    samtools sort "$temp_bam" -o "$out_bam"
    samtools index "$out_bam"
    rm "$temp_bam"
done

