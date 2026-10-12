# Directory containing trimmed .gz files
input_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/RawProcessing"
output_dir="$input_dir/bams"
mkdir -p "$output_dir"

# Reference genome
reference_genome="/mnt/d/xQTL_2025_Data/ref/dm6.fa"

# Number of threads for BWA
threads=50

# Iterate over all R1 files in the specified directory
for file1 in "$input_dir"/*_merged_R1_*_val_1.fq.gz; do
    file2="${file1/_R1_/_R2_}"
    file2="${file2/_val_1/_val_2}"
    echo "$file1"
    echo "$file2"

    if [ -f "$file2" ]; then
        output_file1=$(basename "${file1%_R1_*.fq.gz}_temp_aligned.bam")
        output_file2=$(basename "${file1%_R1_*.fq.gz}_aligned.bam")
        bwa mem "$reference_genome" -M -t "$threads" -v 3 "$file1" "$file2" | samtools view -bS - > "$output_dir/$output_file1"
        samtools sort "$output_dir/$output_file1" -o "$output_dir/$output_file2"
        samtools index "$output_dir/$output_file2"
    else
        echo "Corresponding R2 file for $file1 not found."
    fi
done

