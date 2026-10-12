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
  
  f_mat <- as.matrix(founder_cols[complete.cases(founder_cols), ])
  if (nrow(f_mat) < 10 | ncol(f_mat) < 2) {
    return(tibble(
      sample = pool_names,
      Groups = list(NA), Haps = list(NA), 
      Err = list(NA), Names = list(NA)
    ))
  }
  
  d_mat      <- dist(t(f_mat))
  h_adaptive <- quantile(as.vector(d_mat), 0.90)
  Groups     <- cutree(hclust(d_mat), h = h_adaptive)
  cat(sprintf("SNPs after complete.cases: %d, n_groups: %d, h_p90: %.2f, max_dist: %.2f\n", 
              nrow(f_mat), length(unique(Groups)), h_adaptive, max(d_mat)))
  
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
  founder_mat <- sampdf %>% select(-freq)
  Y           <- sampdf$freq
  
  good        <- !is.na(Y) & complete.cases(founder_mat)
  Y           <- Y[good]
  founder_mat <- as.matrix(founder_mat[good, ])
  
  if (nrow(founder_mat) < 10 | ncol(founder_mat) < 2) {
    return(list(Groups = NA, Haps = NA, Err = NA, Names = NA))
  }
  
  unique_groups <- sort(unique(Groups))
  n_groups      <- length(unique_groups)
  group_names   <- paste0("G", unique_groups)
  
  # Collapse founders to group means BEFORE solving
  A_collapsed <- do.call(cbind, lapply(unique_groups, function(g) {
    rowMeans(founder_mat[, which(Groups == g), drop = FALSE])
  }))
  colnames(A_collapsed) <- group_names
  
  out <- lsei(
    A = A_collapsed, B = Y,
    E = t(matrix(rep(1, n_groups))), F = 1,
    G = diag(n_groups), H = rep(0, n_groups),
    verbose = FALSE, fulloutput = TRUE
  )
  
  names(out$X) <- group_names
  
  list(
    Groups = Groups,
    Haps   = out$X,
    Err    = out$cov,
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

