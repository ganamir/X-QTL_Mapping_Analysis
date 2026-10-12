library(tidyverse)
library(limSolve)
library(data.table)
library(dtplyr)

mychr <- "chr2L"
myparfile <- "haplotype.parameters.R"
mydir <- "process/"

source(myparfile)

filein  <- paste0(mydir, "/RefAlt.",     mychr, ".txt")
rdsfile <- paste0(mydir, "/R.haps.",     mychr, ".rds")
fileout <- paste0(mydir, "/R.haps.",     mychr, ".out.rds")

cat("=== Chromosome:", mychr, "===\n")
cat("Input file:", filein, "\n")

# =============================================================================
# FUNCTIONS
# =============================================================================

# Estimate haplotype composition for all samples in a window
# Works on pre-pivoted wide-format data (no per-window pivot_wider needed)
est_hap <- function(spotsdf, df_wide) {
  chunk <- df_wide %>%
    filter(CHROM == spotsdf$CHROM &
             POS > spotsdf$start &
             POS < spotsdf$end)
  
  founder_cols <- chunk %>% select(all_of(founders))
  pool_names   <- intersect(names_in_bam, names(chunk))
  
  # Cluster once using founder data only (shared across all samples)
  f_mat <- as.matrix(founder_cols[complete.cases(founder_cols), ])
  if (nrow(f_mat) < 10 | ncol(f_mat) < 2) {
    return(tibble(
      sample = pool_names,
      Groups = list(NA), Haps = list(NA), 
      Err = list(NA), Names = list(NA)
    ))
  }
  Groups <- cutree(hclust(dist(t(f_mat))), h = h_cutoff)
  cat(sprintf("SNPs after complete.cases: %d, n_groups: %d, max_dist: %.2f\n", 
              nrow(f_mat), length(unique(Groups)), max(dist(t(f_mat)))))
  
  results <- lapply(pool_names, function(samp) {
    sampdf <- data.frame(freq = chunk[[samp]], founder_cols)
    est_hap2(sampdf, Groups)
  })
  
  tibble(
    sample = pool_names,
    Groups = lapply(results, `[[`, "Groups"),
    Haps   = lapply(results, `[[`, "Haps"),
    Err    = lapply(results, `[[`, "Err"),
    Names  = lapply(results, `[[`, "Names")
  )
}

# Estimate haplotype proportions for a single sample via lsei
# sampdf has columns: freq (pool), then one column per founder
est_hap2 <- function(sampdf, Groups) {
  founder_mat   <- sampdf %>% select(-freq)
  Y             <- sampdf$freq
  
  good          <- !is.na(Y) & complete.cases(founder_mat)
  Y             <- Y[good]
  founder_mat   <- founder_mat[good, ]
  m_founder_mat <- as.matrix(founder_mat)
  
  if (nrow(m_founder_mat) < 10 | ncol(m_founder_mat) < 2) {
    return(list(Groups = NA, Haps = NA, Err = NA, Names = NA))
  }
  
  d <- ncol(m_founder_mat)
  
  out <- lsei(
    A = m_founder_mat, B = Y,
    E = t(matrix(rep(1, d))), F = 1,
    G = diag(d), H = matrix(rep(0.0003, d)),
    verbose = FALSE, fulloutput = TRUE
  )
  
  # Collapse using the shared clustering
  unique_groups <- sort(unique(Groups))
  n_groups      <- length(unique_groups)
  
  T_mat <- matrix(0, nrow = n_groups, ncol = d)
  for (i in seq_along(unique_groups)) {
    T_mat[i, which(Groups == unique_groups[i])] <- 1
  }
  
  collapsed_haps <- as.vector(T_mat %*% out$X)
  collapsed_err  <- T_mat %*% out$cov %*% t(T_mat)
  group_names    <- paste0("G", unique_groups)
  names(collapsed_haps) <- group_names
  
  list(
    Groups = Groups,
    Haps   = collapsed_haps,
    Err    = collapsed_err,
    Names  = group_names
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("Loading RefAlt file...\n")

df <- lazy_dt(read.table(filein, header = TRUE))

df2 <- df %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(RefAlt = str_sub(lab, 1, 3),
         name   = str_sub(lab, 5)) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(freq = REF / (REF + ALT),
         N    = REF + ALT) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

rm(df)
cat("df2 loaded:", nrow(df2), "rows\n")

# =============================================================================
# FOUNDER FREQUENCY ROUNDING
# =============================================================================

df3 <- df2 %>%
  mutate(freq = case_when(
    name %in% founders & freq >  0.55 ~  1,
    name %in% founders & freq <  0.45 ~  0,
    name %in% founders & freq >= 0.45 & freq <= 0.55 ~ NA_real_,
    TRUE ~ freq
  )) %>%
  filter(!(name %in% founders & is.na(freq)))

rm(df2)
cat("df3 after founder rounding:", nrow(df3), "rows\n")

# =============================================================================
# SAVE PROCESSED DATA (long format, for downstream use)
# =============================================================================

saveRDS(df3, file = rdsfile)

# =============================================================================
# PIVOT TO WIDE FORMAT ONCE
# =============================================================================

cat("Pivoting to wide format (one-time operation)...\n")

df3_wide <- df3 %>%
  select(-N) %>%
  pivot_wider(names_from = name, values_from = freq)

cat("df3_wide:", nrow(df3_wide), "rows x", ncol(df3_wide), "cols\n")

# --- DIAGNOSTIC ---
test_chunk <- df3_wide %>% 
  filter(CHROM == mychr) %>% 
  slice(1:500)

founder_cols <- test_chunk %>% select(all_of(founders))
cat("Total SNPs in chunk:", nrow(founder_cols), "\n")
cat("SNPs after complete.cases:", sum(complete.cases(founder_cols)), "\n")

f_mat <- as.matrix(founder_cols[complete.cases(founder_cols), ])
if (nrow(f_mat) > 1) {
  d <- max(dist(t(f_mat)))
  cat("Max pairwise distance:", d, "\n")
  cat("Groups at h=7.5:", length(unique(cutree(hclust(dist(t(f_mat))), h = 7.5))), "\n")
}


# Test with an actual window, not just first 500 rows
mid_pos <- median(df3_wide$POS[df3_wide$CHROM == mychr])
test_chunk <- df3_wide %>% 
  filter(CHROM == mychr & POS > (mid_pos - size) & POS < (mid_pos + size))

founder_cols <- test_chunk %>% select(all_of(founders))
cat("Total SNPs in window:", nrow(founder_cols), "\n")
cat("SNPs after complete.cases:", sum(complete.cases(founder_cols)), "\n")

f_mat <- as.matrix(founder_cols[complete.cases(founder_cols), ])
d <- max(dist(t(f_mat)))
cat("Max pairwise distance:", d, "\n")

# Check the distribution of distances
all_d <- as.vector(dist(t(f_mat)))
cat("Distance quantiles:\n")
print(quantile(all_d, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 1.0)))

# Try a few cutoffs
for (h in c(0.5, 1.0, 1.5, 2.0)) {
  g <- length(unique(cutree(hclust(dist(t(f_mat))), h = h)))
  cat(sprintf("h=%.1f: %d groups\n", h, g))
}

# --- END DIAGNOSTIC ---

rm(df3)
gc()

# =============================================================================
# DEFINE SCAN WINDOWS
# =============================================================================

minpos <- min(df3_wide$POS)
maxpos <- max(df3_wide$POS)
myseq  <- seq(0, maxpos, step)
myseq  <- myseq[myseq > minpos + size & myseq < maxpos - size]

spots <- data.frame(
  CHROM = rep(mychr, length(myseq)),
  pos   = myseq,
  start = myseq - size,
  end   = myseq + size
)

UU    <- unique(df3_wide$POS)
spots <- spots %>%
  rowwise() %>%
  mutate(NN = sum(start < UU) - sum(end < UU)) %>%
  filter(NN >= 50) %>%
  select(-NN)

cat(sprintf("Scan windows: %d\n", nrow(spots)))

# =============================================================================
# HAPLOTYPE SCAN (sequential)
# =============================================================================

total_windows <- nrow(spots)
cat(sprintf("Starting scan: %d windows, h_cutoff = %.1f\n", total_windows, h_cutoff))

t_start <- Sys.time()

spots2 <- spots %>%
  group_nest(row_number()) %>%
  mutate(out = imap(data, function(x, i) {
    if (i %% 100 == 0) {
      cat(sprintf("Progress: %d/%d windows (%.1f%%)\n",
                  i, total_windows, 100 * i / total_windows))
      flush.console()
    }
    est_hap(x, df3_wide)
  })) %>%
  unnest(data) %>%
  select(-c(start, end)) %>%
  select(-`row_number()`) %>%
  unnest_wider(out)

t_end <- Sys.time()
cat(sprintf("Scan finished in %.1f minutes\n", difftime(t_end, t_start, units = "mins")))

saveRDS(spots2, file = fileout)
cat("Output saved:", fileout, "\n")





