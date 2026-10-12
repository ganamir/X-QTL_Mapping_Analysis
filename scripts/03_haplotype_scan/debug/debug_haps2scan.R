# =============================================================================
# DEBUG HAPLOTYPEINF
# =============================================================================

design.df <- read.table("input_table.txt")

#########
# Functions
#########

average_variance <- function(cov_matrix, tolerance = 1e-10) {
  n <- nrow(cov_matrix)  
  # Calculate eigenvalues
  eigenvalues <- eigen(cov_matrix, only.values = TRUE)$values  
  # Filter out eigenvalues that are effectively zero or negative
  positive_eigenvalues <- eigenvalues[eigenvalues > tolerance]  
  # Calculate the product of positive eigenvalues
  log_det <- sum(log(positive_eigenvalues))  
  # Use the number of positive eigenvalues for the root
  n_positive <- length(positive_eigenvalues)  
  # Calculate log of average variance
  log_avg_var <- log_det / n_positive  
  # Convert back to original scale
  avg_var <- exp(log_avg_var)  
  return(list(avg_var = avg_var, n_positive = n_positive, n_total = n))
}

wald.test3 = function(p1,p2,covar1,covar2,nrepl=1,N1=NA,N2=NA){
  
  # Wald test for multinomial frequencies
  # if nrepl = 1: (one replicate, analogous to chi square):
  #  p1 and p2 are vectors of relative frequencies to be compared
  # covar1 and covar2 are the reconstruction error 
  # covariance matrices from limSolve
  # the sampling covariance matrices are generated within limSolve
  # if nrepl > 1 (multiple replicates, analogous to CMH):
  #   p1 and p2 are matrices, each row is frequency vector for one replicate
  # covar1 and covar2 are tensors (3-dimensional arrays, third dimension 
  #  denotes replicate) for the linSolve covariance matrices
  # N1 (initial) and N2 (after treatment) 
  # are sample sizes, they are vectors when there is more than one replicate
  # N1[i], N2[i] are then for replicate i
  if (nrepl>1){
    N1.eff=rep(NA,nrepl)
    N2.eff=rep(NA,nrepl)
    lp1 = length(p1[1,])
    cv1=array(NA,c(lp1,lp1,nrepl))
    cv2=array(NA,c(lp1,lp1,nrepl))
    for (i in 1:nrepl){
      
      covmat1  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N1[i])
      covmat2  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N2[i])
      
      N1.eff[i] = sum(diag(covmat1))*4*N1[i]^2/(sum(diag(covmat1))*2*N1[i]+2*N1[i]*sum(diag(covar1[,,i])) )
      N2.eff[i] = sum(diag(covmat2))*4*N2[i]^2/(sum(diag(covmat2))*2*N2[i]+2*N2[i]*sum(diag(covar2[,,i])) )
      cv1[,,i]= (covmat1 + covar1[ , ,i])  * (N1.eff[i])^2
      cv2[,,i]= (covmat2 + covar2[ , ,i])  * (N2.eff[i])^2
      
    }
    
    p1 = N1.eff %*% p1 / sum(N1.eff)
    p2 = N2.eff %*% p2 / sum(N2.eff)
    
    covar1= rowSums(cv1, dims = 2) / sum(N1.eff)^2
    covar2= rowSums(cv2, dims = 2) / sum(N2.eff)^2
    # browser()
  }
  else {
    covmat1  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N1)
    covmat2  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N2)
    covar1 = covar1 + covmat1
    covar2 = covar2 + covmat2
  }
  
  df = length(p1)-1
  covar=covar1+covar2
  eg <- eigen(covar)
  # remove last eigenvector which corresponds to eigenvalue zero
  keep <- eg$values[1:df] > 1e-4  # threshold for meaningful eigenvalues
  ev <- eg$vectors[,1:df][,keep]
  eval <- eg$values[1:df][keep]
  df <- sum(keep)  # adjust degrees of freedom
  trafo<-diag(1/sqrt(eval)) %*% t(ev) 
  # set extremely small values to zero
  #new.covar[new.covar < 10^-9]=0
  p1= as.vector(p1); p2=as.vector(p2)
  tstat <- sum((trafo %*% (p1 - p2))^2)
  pval<- exp(pchisq(tstat,df,lower.tail=FALSE,log.p=TRUE))
  list(wald.test=tstat, p.value=pval, avg.var=average_variance(covar)$avg_var)
}

mn.covmat= function(p,n,min.p=0){
  # generate multinomial covariance matrix
  # p is vector of multinomial relative frequencies
  # n is sample size
  # compute covariance matrix for relative frequencies, for absolute frequencies multiply by n^2
  # if min.p >0, then values of p smaller than min.p are set to min.p and the resulting vector is rescaled.
  p[p<min.p] = min.p; p=p/sum(p)
  mat = - tcrossprod(p)
  diag(mat) = p*(1-p)
  mat = mat/n
  mat
}

pseudoN.test = function(p1,p2,covar1,covar2,nrepl,N1,N2){
  pseudoN_C = rep(NA,nrepl)
  pseudoN_Z = rep(NA,nrepl)
  for(i in 1:nrepl){
    pseudoN_C[i] = (2 * N1[i] * sum(p1[i,] * (1-p1[i]))) / (2 * N1[i] * sum(diag(covar1[,,i])) + sum(p1[i,] * (1-p1[i])))
    pseudoN_Z[i] = (2 * N2[i] * sum(p2[i,] * (1-p2[i]))) / (2 * N2[i] * sum(diag(covar2[,,i])) + sum(p2[i,] * (1-p2[i])))
  }
  Count1 = round(p1*pseudoN_C,0)
  Count2 = round(p2*pseudoN_Z,0)
  lowCountFounder = apply(rbind(Count1,Count2),2,sum)
  if(sum(lowCountFounder>=5)<2){
    log10p = NA
  }else{
    Count1 = Count1[,lowCountFounder >= 5]		
    Count2 = Count2[,lowCountFounder >= 5]		
    if(nrepl==1){
      out=chisq.test(rbind(Count1,Count2),correct=TRUE)
    }else{
      nF = ncol(Count1)
      tdf = data.frame(Count=c(as.numeric(t(Count1)),as.numeric(t(Count2))),
                       founder=rep(1:nF,2*nrepl),
                       TRT = c(rep(1,nF*nrepl),rep(2,nF*nrepl)),
                       REP = c(rep(1:nrepl,each=nF),rep(1:nrepl,each=nF)))
      D.x = xtabs(Count ~ founder + TRT + REP, data = tdf)
      out = mantelhaen.test(D.x,correct=TRUE)
    }
    log10p = -log10(out$p.value)
  }
  log10p
}

add_genetic = function(df){
  df$cM = rep(NA,nrow(df))
  fm=read.table("flymap.r6.txt",header=FALSE)
  colnames(fm)=c("chr","pos","cM")
  library(splines)
  for(chrs in c("chrX","chr2L","chr2R","chr3L","chr3R")){
    fmX = fm %>% filter(chr==chrs)
    out = ksmooth(fmX$pos,fmX$cM,kernel="normal",bandwidth=3e6)
    f_of_x = splinefun(out$x,out$y)
    temp = f_of_x(df$pos[df$chr==chrs])
    df$cM[df$chr==chrs] = temp
  }
  df
}

Heritability = function(p1, p2, nrepl, ProportionSelect, af_cutoff){
  nF = ncol(p1)
  tdf = data.frame(freq=c(as.numeric(t(p1)),as.numeric(t(p2))),
                   founder=rep(1:nF,2*nrepl),
                   TRT = c(rep("C",nF*nrepl),rep("Z",nF*nrepl)),
                   REP = c(rep(1:nrepl,each=nF),rep(1:nrepl,each=nF)))
  
  Falconer_H2 = tdf %>%
    pivot_wider(names_from = TRT, values_from = freq) %>%
    mutate(mean_diff_sq = (Z-C)^2) %>%
    mutate(mean_af_C = case_when(C <= af_cutoff ~ af_cutoff, .default = C)) %>%
    mutate(H2temp = mean_diff_sq/mean_af_C) %>%
    group_by(REP) %>%
    summarize(H2temp_sum = sum(H2temp)) %>%
    ungroup() %>%
    left_join(ProportionSelect,by="REP") %>%
    filter(!is.na(Proportion)) %>%
    mutate(Falcon_i = dnorm(qnorm(1-Proportion))/Proportion) %>%
    group_by(REP) %>%
    summarize(H2 = 200 * H2temp_sum / Falcon_i^2) %>%
    ungroup() %>%
    summarize(mH2 = mean(H2)) %>%
    pull(mH2)
  
  Cutler_H2 = tdf %>%
    pivot_wider(names_from = TRT, values_from = freq) %>%
    left_join(ProportionSelect,by="REP") %>%
    filter(!is.na(Proportion)) %>%
    mutate(Penetrance = (Z * Proportion)/C) %>%
    mutate(Penetrance = case_when(Penetrance <= Proportion/2 ~ Proportion/2,
                                  Penetrance >= 2*Proportion ~ 2*Proportion,
                                  .default = Penetrance)) %>% 
    mutate(Affect = qnorm(1-Proportion) - qnorm(1-Penetrance)) %>%
    mutate(marg_Va = Affect^2 * C) %>%
    group_by(REP) %>%
    mutate(H2 = 200*sum(marg_Va)) %>%
    ungroup() %>%
    summarize(mH2 = mean(H2)) %>%
    pull(mH2)
  
  list(Falconer_H2=Falconer_H2, Cutler_H2=Cutler_H2)
}

doscan = function(df,chr,Nfounders){
  sexlink = 1
  if(chr=="chrX"){ sexlink=0.75 }
  
  # I tested with xx2$data[[1]]
  df2 = df %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(recodeTable) %>%
    select(-sample) %>% mutate(sample=pool) %>% select(-pool) %>%
    left_join(Numflies, join_by(sample==pool)) %>%
    separate(sample,into=c("longTRT","REP","REPrep"),remove=FALSE) %>%
    left_join(TreatmentMapping)
  
  # only analyze data for which all founders are discernable..
  allFounders = as.numeric(df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm)))	
  
  ll = list(Wald_log10p = NA, Pseu_log10p = NA, Falc_H2 = NA, Cutl_H2 = NA, avg.var = NA)
  if(allFounders!=Nfounders){ return(ll) }
  
  ##  now cases where all founders are OK
  ##  now collapse any pure replicates.  This is tidy ugly.  But I feel there 
  ##  is value in keeping dataframe columns as lists...
  df3 = df2 %>%
    select(-Groups) %>%
    group_by(TRT,REP) %>%	
    summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
              Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
              Names = list(first(Names)),
              Num_mean = sexlink*mean(Num)) %>%
    rename(Haps=Haps_mean,Num=Num_mean,Err=Err_mean)
  
  ## these summaries of the data are pretty useful for tests
  p1 = df3 %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t() 
  row.names(p1) <- NULL
  p2 = df3 %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  row.names(p2) <- NULL
  covar1 = do.call(abind, c(df3 %>% filter(TRT=="C") %>% pull(Err), along = 3))
  covar2 = do.call(abind, c(df3 %>% filter(TRT=="Z") %>% pull(Err), along = 3))
  nrepl = df3 %>% filter(TRT=="C") %>% nrow()
  nrepl == df3 %>% filter(TRT=="Z") %>% nrow()
  N1 = df3 %>% filter(TRT=="C") %>% pull(Num)
  N2 = df3 %>% filter(TRT=="Z") %>% pull(Num)
  
  wt = tryCatch(
    wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) list(wald.test=NA, p.value=NA, avg.var=NA)
  )
  Wald_log10p = -log10(wt$p.value)
  Pseu_log10p = tryCatch(
    pseudoN.test(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) NA
  )
  
  af_cutoff = 0.01     # 1% --- heritability estimators can be off for really low allele frequencies
  temp = Heritability(p1, p2, nrepl, ProportionSelect, af_cutoff)
  Falc_H2 = temp$Falconer_H2
  Cutl_H2 = temp$Cutler_H2
  
  ll = list(Wald_log10p = Wald_log10p, Pseu_log10p = Pseu_log10p,
            Falc_H2 = Falc_H2, Cutl_H2 = Cutl_H2, avg.var = wt$avg.var)
  ll
}

doscan2 = function(df,chr,Nfounders){
  sexlink = 1
  if(chr=="chrX"){ sexlink=0.75 }
  
  # I tested with xx2$data[[1]]
  df2 = df %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(design.df, join_by(sample==bam)) %>%
    filter(!is.na(TRT))
  
  # only analyze data for which all founders are discernable..
  allFounders = as.numeric(df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm)))	
  
  ll = list(Wald_log10p = NA, Pseu_log10p = NA, Falc_H2 = NA, Cutl_H2 = NA, avg.var = NA)
  if (is.na(allFounders) || allFounders < 2) { return(ll) }
  
  ##  now cases where all founders are OK
  ##  now collapse any pure replicates.  This is tidy ugly.  But I feel there 
  ##  is value in keeping dataframe columns as lists...
  df3 = df2 %>%
    select(-Groups) %>%
    group_by(TRT,REP) %>%	
    summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
              Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
              Names = list(first(Names)),
              Num_mean = sexlink*mean(Num)) %>%
    rename(Haps=Haps_mean,Num=Num_mean,Err=Err_mean)
  
  ## these summaries of the data are pretty useful for tests
  p1 = df3 %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t() 
  row.names(p1) <- NULL
  p2 = df3 %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  row.names(p2) <- NULL
  covar1 = do.call(abind, c(df3 %>% filter(TRT=="C") %>% pull(Err), along = 3))
  covar2 = do.call(abind, c(df3 %>% filter(TRT=="Z") %>% pull(Err), along = 3))
  nrepl = df3 %>% filter(TRT=="C") %>% nrow()
  nrepl == df3 %>% filter(TRT=="Z") %>% nrow()
  N1 = df3 %>% filter(TRT=="C") %>% pull(Num)
  N2 = df3 %>% filter(TRT=="Z") %>% pull(Num)
  
  wt = tryCatch(
    wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) list(wald.test=NA, p.value=NA, avg.var=NA)
  )
  Wald_log10p = -log10(wt$p.value)
  Pseu_log10p = tryCatch(
    pseudoN.test(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) NA
  )
  
  af_cutoff = 0.01     # 1% --- heritability estimators can be off for really low allele frequencies
  temp = Heritability(p1, p2, nrepl, ProportionSelect, af_cutoff)
  Falc_H2 = temp$Falconer_H2
  Cutl_H2 = temp$Cutler_H2
  
  ll = list(Wald_log10p = Wald_log10p, Pseu_log10p = Pseu_log10p,
            Falc_H2 = Falc_H2, Cutl_H2 = Cutl_H2, avg.var = wt$avg.var)
  ll
}

xx1 = readRDS("process/R.haps.chr2L.out.rds")
Nfounders=length(xx1$Groups[[1]][[1]])
ProportionSelect = design.df %>% filter(TRT=="Z") %>% select(REP,Proportion) %>% arrange(REP)

bb1 = xx1 %>%
  head(n=100) %>%
  group_by(CHROM,pos) %>%
  nest() %>%
  mutate(out = map2(data, CHROM, doscan2, Nfounders=Nfounders)) %>%
  unnest_wider(out)
bb2 = bb1 %>% select(-data) %>% rename(chr=CHROM)
bb3 = add_genetic(bb2)

# ---- DIAGNOSTIC: inspect eigenvalues on a bad window ----
if (any(bb3$Wald_log10p > 200, na.rm = TRUE)) {
  bad_row <- bb3 %>% filter(Wald_log10p > 200) %>% slice(1)
  bad_pos <- bad_row$pos
  bad_chr <- bad_row$chr
  cat("Inspecting bad window at", bad_chr, "pos =", bad_pos, "\n")
  
  bad_data <- xx1 %>% filter(pos == bad_pos)
  
  df2_bad <- bad_data %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(design.df, join_by(sample==bam)) %>%
    filter(!is.na(TRT))
  
  df3_bad <- df2_bad %>%
    select(-Groups) %>%
    group_by(TRT,REP) %>%
    summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
              Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
              Names = list(first(Names)),
              Num_mean = mean(Num)) %>%
    rename(Haps=Haps_mean, Num=Num_mean, Err=Err_mean)
  
  p1_bad = df3_bad %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  p2_bad = df3_bad %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  N1_bad = df3_bad %>% filter(TRT=="C") %>% pull(Num)
  N2_bad = df3_bad %>% filter(TRT=="Z") %>% pull(Num)
  covar1_bad = do.call(abind, c(df3_bad %>% filter(TRT=="C") %>% pull(Err), along = 3))
  covar2_bad = do.call(abind, c(df3_bad %>% filter(TRT=="Z") %>% pull(Err), along = 3))
  nrepl_bad = df3_bad %>% filter(TRT=="C") %>% nrow()
  
  N1.eff = N2.eff = rep(NA, nrepl_bad)
  lp1 = ncol(p1_bad)
  cv1 = cv2 = array(NA, c(lp1, lp1, nrepl_bad))
  for (i in 1:nrepl_bad) {
    covmat1 = mn.covmat((N1_bad[i]*p1_bad[i,]+N2_bad[i]*p2_bad[i,])/(N1_bad[i]+N2_bad[i]), 2*N1_bad[i])
    covmat2 = mn.covmat((N1_bad[i]*p1_bad[i,]+N2_bad[i]*p2_bad[i,])/(N1_bad[i]+N2_bad[i]), 2*N2_bad[i])
    N1.eff[i] = sum(diag(covmat1))*4*N1_bad[i]^2/(sum(diag(covmat1))*2*N1_bad[i]+2*N1_bad[i]*sum(diag(covar1_bad[,,i])))
    N2.eff[i] = sum(diag(covmat2))*4*N2_bad[i]^2/(sum(diag(covmat2))*2*N2_bad[i]+2*N2_bad[i]*sum(diag(covar2_bad[,,i])))
    cv1[,,i] = (covmat1 + covar1_bad[,,i]) * (N1.eff[i])^2
    cv2[,,i] = (covmat2 + covar2_bad[,,i]) * (N2.eff[i])^2
  }
  covar_total = rowSums(cv1, dims = 2)/sum(N1.eff)^2 + rowSums(cv2, dims = 2)/sum(N2.eff)^2
  
  evals = eigen(covar_total)$values
  cat("Number of clusters:", lp1, "\n")
  cat("All eigenvalues:\n")
  print(signif(sort(evals, decreasing = TRUE), 4))
  cat("\nEigenvalue range:", signif(max(evals), 4), "to", signif(min(evals), 4), "\n")
  stop("Diagnostic done - remove this block and set proper threshold")
} else {
  cat("No windows with Wald_log10p > 200 found\n")
}

# Add to the diagnostic, after covar_total:
cat("\nDiag of lsei covariance (reconstruction error):\n")
print(signif(diag(covar1_bad[,,1]), 4))
cat("\nDiag of multinomial covariance (sampling error):\n")
pooled_p <- (N1_bad[1]*p1_bad[1,]+N2_bad[1]*p2_bad[1,])/(N1_bad[1]+N2_bad[1])
print(signif(diag(mn.covmat(pooled_p, 2*N1_bad[1])), 4))



cat("\nHaps for Control rep 1:\n")
print(signif(p1_bad[1,], 4))
cat("\nHaps for Treatment rep 1:\n")
print(signif(p2_bad[1,], 4))
cat("\nDifference:\n")
print(signif(p2_bad[1,] - p1_bad[1,], 4))



# Check how many windows have out-of-range proportions
bad_count <- 0
for (i in 1:nrow(xx1)) {
  haps <- xx1$Haps[[i]]
  if (!is.list(haps)) next
  for (h in haps) {
    if (any(!is.na(h) & (h < -0.01 | h > 1.01))) {
      bad_count <- bad_count + 1
      break
    }
  }
}
cat("Windows with invalid proportions:", bad_count, "out of", nrow(xx1), "\n")



# I drop the loci for which the scan gives and NA
#bb4 = bb1 %>%
#	filter(!is.na(Pseu_log10p)) %>%
#	select(-c(Wald_log10p, Pseu_log10p, Falc_H2, Cutl_H2, avg.var, data)) %>%
#	left_join(xx1) %>%
#	select(-c(Err,Groups)) %>%
#	unnest(c(sample,Haps,Names)) %>%
#	unnest(c(Haps,Names)) %>%
#	rename(chr=CHROM,pool=sample,freq=Haps,founder=Names) %>%
#	left_join(design.df, by=c("pool"="bam")) %>%
#	select(c(chr,pos,TRT,REP,REPrep,freq,founder)) %>%
#	filter(!is.na(TRT)) %>%
#	group_by(chr,pos,TRT,REP,founder) %>%
#	summarize(freq=mean(freq,na.rm=TRUE))

write.table(bb3, fileout)
#write.table(bb4, fileout_meansBySample)
