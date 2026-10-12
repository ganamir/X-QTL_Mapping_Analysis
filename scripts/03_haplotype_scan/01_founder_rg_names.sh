for bam in /mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams_RG/DGRP_*.bam; do
    samtools view -H "$bam" | grep '^@RG' | grep -oP 'SM:\S+' | sed 's/SM://'
done | sort -V | awk 'BEGIN{printf "c("} {printf sep"\""$0"\""; sep=","} END{print ")"}' \
> /mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/founder_RG.txt
