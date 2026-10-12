# X-QTL Mapping Analysis

Code for **"Repeated genotype–phenotype mapping reveals context dependency and the limits of genetic maps for predicting genomic adaptation"** (Gabidulin et al.).

Preprint: *bioRxiv DOI to be added* · Archived release: *Zenodo DOI to be added*

## Overview

Spinosad resistance in *Drosophila melanogaster* was mapped four times by extreme-QTL (X-QTL), changing one factor at a time, and compared with genomic change in a field evolve-and-resequence (E&R) population exposed to spinosad.

| Experiment | Mapping population | Dose | Environment |
|---|---|---|---|
| LC50-DSPR-I (reference) | DSPR | LC50 | Indoor |
| LC50-DGRP-I | DGRP | LC50 | Indoor |
| LC90-DSPR-I | DSPR | LC90 | Indoor |
| LC50-DSPR-O | DSPR | LC50 | Outdoor |
| Field E&R | Outdoor mesocosms | Exposed vs control | Outdoor, 4 collections (Aug 24 to Nov 24) |

## Repository layout

```
scripts/
  01_dgrp_founders/     download, trim, align, dedup and add read groups for the ~100 DGRP founder lines
  02_xqtl_pools/        process X-QTL pool sequencing (selected and control replicates)
  03_haplotype_scan/    founder haplotype inference and X-QTL scans (modified XQTL2 pipeline)
    modified_DGRP/      XQTL2 scripts modified to handle ~100 DGRP founders
    debug/              diagnostics used to choose clustering and scan settings
notebooks/              downstream analysis (R, Jupyter), run in order 01 to 04, saved with outputs
```

Input data and intermediate files are not stored here; they are archived on Zenodo (see Data).

Scripts contain the absolute paths used on our system (`/mnt/d/...`). Edit the directory variables at the top of each script before running.

## Pipeline

### 1. DGRP founder preparation (`scripts/01_dgrp_founders/`)

1. Download DGRP line information from [DGRPool](https://dgrpool.epfl.ch/) and filter to the founder lines used in the DGRP mapping population.
2. Match lines to the Huang et al. (2014) resequencing data ([supplement](https://genome.cshlp.org/content/24/7/1193/suppl/DC1)), keeping Illumina runs only (454 runs dropped).
3. `01_download_sra.sh`: download each line's SRA runs and convert to FASTQ (SLURM array).
4. `02_trim.sh`: trim with Trim Galore (`-q 20 --length 30 --max_n 1`), handling paired- and single-end runs.
5. `03_align.sh`: align to dm6 with `bwa mem`, then sort and index.
6. `04_dedup_merge.sh`: remove duplicates per run (`samtools markdup -r`), then merge runs per line.
7. `05_add_readgroups.sh`: set one read group per line (`DGRP_<n>`).

### 2. X-QTL pool processing (`scripts/02_xqtl_pools/`)

1. `01_merge_lanes_DGRP_only.sh`: merge lanes (DGRP X-QTL only; DSPR starts at step 2).
2. `02_trim.sh`: Trim Galore (`-q 20 --length 40 --max_n 1`).
3. `03_align.sh`: `bwa mem` to dm6, then sort and index.
4. `04_dedup.sh`: duplicate removal with `samtools markdup -r`.
5. `05_add_readgroups.sh`: rename pools to `A1`–`A8` (selected) and `C1`–`C8` (control), add read groups and index.

### 3. Haplotype inference and scans (`scripts/03_haplotype_scan/`)

Built on the Long lab's [XQTL2](https://github.com/tdlong/XQTL2) pipeline.

1. `01_founder_rg_names.sh`: collect founder read-group names for `haplotype.parameters.R`.
2. `02_refalt_tables.sh`: joint `bcftools mpileup`/`call` on founders plus pools, per chromosome; biallelic SNPs with QUAL > 59 written as REF/ALT count tables.
3. Edit `haplotype.parameters.R` (founders, `step`, `size`, `h_cutoff`, pool names) and `input_table.txt` (replicate, treatment, pool size and proportion selected).
4. Run `REFALT2haps.Andreas.sh haplotype.parameters.R process/`: estimates founder haplotype frequencies in each pool in windows of ±`size` around scan positions spaced `step` apart, by least-squares estimation under inequality constraints (LSEI; frequencies non-negative and summing to one).
5. Run `haps2scan.Apr2025.sh input_table.txt process/ <name>`: Wald test of haplotype frequency differences between selected and control pools per scan position.

**Modifications for the DGRP mapping population** (`modified_DGRP/`). Founder coverage is low (~5–50×) and there are ~100 founders, so the original scripts dropped most sites and produced unstable fits. Changes:

- **`REFALT2haps.Andreas.code.r`**
  - Founder allele frequencies are binarized (> 0.55 → 1, < 0.45 → 0, 0.45–0.55 → missing), reflecting the expectation that inbred lines are homozygous.
  - Founders are clustered in each window (`hclust`, cut at the 90th percentile of pairwise distances, giving about 6–10 groups, similar to the 8 DSPR founders), and LSEI is solved on group means rather than individual founders. DSPR uses a fixed `h_cutoff = 2.5`.
- **`haps2scan.Apr2025.code.r`**: the scan accepts a variable number of founder groups per window and handles missing values from failed fits.
- **`scan_functions.r`**
  - Covariance eigenvalues below 1e-5 are dropped before the Wald test, as a guard for windows where clustering does not resolve founders.
  - Per-founder contribution output is disabled, since founders are grouped.

### 4. Downstream analysis (`notebooks/`)

Run in order 01 to 04. Each notebook sets its paths in the first cell (`ROOT`, `IN_DIR`, `OUT_DIR`); point `IN_DIR` at the unpacked Zenodo inputs and `OUT_DIR` at an empty folder, where each notebook writes `Data/<notebook>/` and `Figures/<notebook>/` for the next one to read. The notebooks are saved with the outputs used in the paper.

| Notebook | What it does | Paper outputs |
|---|---|---|
| `01_GLMs.ipynb` | Per-SNP quasibinomial GLMs for field exposed vs control (each collection and all collections combined); converts field SNP and X-QTL window *P* to \|*z*\|; builds 50/100/500 kb window tables | Inputs to 02–04 |
| `02_Peak_Analysis.ipynb` | Outlier windows (top 5% per experiment), triple-concordant windows (≥ 3 of 4 experiments), triple regions (scan footprint, ±50 kb), genes, circular-shift Spearman tests, decile decomposition, FlyCADD candidates | Figures 2–4; Figures S5–S8; Table S1; File S1 |
| `03_Supplemental.ipynb` | Selection intensity, permutation tests of ΔAF in triple regions, DGRP vs DSPR site overlap, per-SNP ΔAF across experiments and field | Figures S1, S3, S4, S11; Table S4 |
| `04_Transcriptomic_overlap.ipynb` | Overlap of peak genes with spinosad-resistance DEGs (Tabuloc et al. 2024): hypergeometric and rank tests | Figures S9–S10; Tables S2–S3 |

**Triple regions.** Each X-QTL value at a scan position summarizes ±50 kb, so gene- and SNP-level analyses use the scan footprint of the triple windows (first scan position − 50 kb to last + 50 kb; peaks ≤ 500 kb apart merged, which gives the 2L Consensus region), not the 50 kb windows. This keeps gene and SNP sets stable when the 50 kb sampling grid shifts.

## Data

- **Raw reads:** NCBI SRA, BioProject *PRJNA to be added*.
- **Notebook inputs and processed outputs:** Zenodo, *DOI to be added*. The archive is sorted by type and experiment, with a `MANIFEST.csv` giving each file's original path, size and md5.
- **DGRP founder reads:** Huang et al. (2014), via SRA.
- **FlyCADD scores:** from the FlyCADD authors; not redistributed here.

## Requirements

- **Upstream:** SRA Toolkit, Trim Galore, BWA, SAMtools, Picard, BCFtools, and SLURM for the array scripts.
- **R ≥ 4.4:** `tidyverse`, `data.table`, `dtplyr`, `limSolve`, `abind`, `ggplot2`, `patchwork`, `ggtext`, `ggrastr`, `scales`, `readxl`, `doSNOW`, `foreach`.

## Citation

If you use this code, please cite the paper (citation to be added) and the original XQTL2 pipeline (Long lab, https://github.com/tdlong/XQTL2).

## License

*To be added.*
