# JUICE Haplotype Parameters
# This file defines the parameters for haplotype estimation in the JUICE dataset

# Founder samples (these should match the column names in your REFALT data)
founders <- c("DGRP_21","DGRP_26","DGRP_31","DGRP_48","DGRP_69","DGRP_73","DGRP_83","DGRP_85","DGRP_88","DGRP_91","DGRP_100","DGRP_101","DGRP_105","DGRP_109","DGRP_129","DGRP_136","DGRP_142","DGRP_149","DGRP_153","DGRP_161","DGRP_176","DGRP_177","DGRP_181","DGRP_189","DGRP_217","DGRP_227","DGRP_228","DGRP_239","DGRP_287","DGRP_301","DGRP_306","DGRP_309","DGRP_310","DGRP_313","DGRP_315","DGRP_318","DGRP_321","DGRP_324","DGRP_332","DGRP_335","DGRP_338","DGRP_340","DGRP_350","DGRP_356","DGRP_359","DGRP_365","DGRP_370","DGRP_371","DGRP_373","DGRP_374","DGRP_375","DGRP_380","DGRP_381","DGRP_382","DGRP_385","DGRP_386","DGRP_390","DGRP_391","DGRP_399","DGRP_405","DGRP_406","DGRP_427","DGRP_437","DGRP_440","DGRP_441","DGRP_443","DGRP_502","DGRP_509","DGRP_513","DGRP_517","DGRP_535","DGRP_551","DGRP_559","DGRP_563","DGRP_566","DGRP_589","DGRP_634","DGRP_703","DGRP_705","DGRP_712","DGRP_732","DGRP_737","DGRP_748","DGRP_757","DGRP_776","DGRP_783","DGRP_786","DGRP_801","DGRP_804","DGRP_812","DGRP_819","DGRP_820","DGRP_821","DGRP_852","DGRP_857","DGRP_859","DGRP_861","DGRP_890","DGRP_900","DGRP_908")
#founders <- readLines("founder_names.txt")
# Step size for scanning positions along chromosome
step <- 1000  # 10kb steps (increased from 5kb for speed)

#try 1kb, 50kb, 100kb, 200kb

#maybe 10kb

# Window size - COMMENTED OUT, comes from command line
size <- 50000  # Window size comes from command line parameter

#try 100kb, 200kb, 300kb, 400kb, 500kb

# Clustering parameters
h_cutoff <- 2.5  # Default h_cutoff for fixed window hierarchical clustering

# Samples to process (allows selecting subset from large REFALT files)
names_in_bam=c("C1","C2","C3","C4","C5","C6","C7","C8","A1","A2","A3","A4","A5","A6","A7","A8")

