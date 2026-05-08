####################################################################
# Who: Amanda Mejia
# When: Fall 2024-2025
# What: Compute effective sample size and baseline variance
# How: First, 4_FC_ACFs.R, 4_AggFCandFlags.R, 5_AggFlags.R, 5_error.R, and 6_ScrubbingRate.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## ---------------------------------------------------------------------------
source("code/0_SharedCode.R")
setwd(dir_slate)

library(dplyr) #version 1.1.4

gt <- "FD2"

A_sessions <- c("test_LR1","test_RL1")
B_sessions <- c("test_LR2","test_RL2")

aggFD <- readRDS(file.path(dir_results, "results/6_scrubRate/AggFD_withSubjectSplit.rds"))
load(file = file.path(dir_github,'data','subjects_retest.RData')) #sub_drop, sub_himo, sub_lomo (from 6_ScrubbingRate.R)
#sub_drop based on ground truth scan duration < 50 min

# Scrubbing Flags --------------------------------------------------------------

# Collect scrubbing rate by method (T=10 min) 
flags <- readRDS(file.path(dir_results, 'results', '5_flags', 'withDVARS_flags.rds'))
saveRDS(flags, file.path(dir_github, 'results', '5_flags', 'withDVARS_flags.rds')) #copy to Github
dimnames(flags)[[5]] <- c("None","FD5","FD4","FD3","FD2")
flags_ <- flags[,,,,,,"FD", drop = FALSE]
flags_plus <- flags[,,,,"FD2",,"FDplus", drop = FALSE]
dimnames(flags_plus)[[5]] <- "FD2+"
flags <- abind::abind(flags_, flags_plus, along = 5)[,,,,,,1]
rm(flags_, flags_plus)

#number of volumes scrubbed at T = 5 min
flags_5min <- flags[c("test_LR1","test_LR2","test_RL1","test_RL2"),,,"first_2.5_mins","FD2+",1]
flags_5min <- apply(flags_5min, 1:2, sum, na.rm = TRUE) #count the total number of volumes scrubbed for each test
flags_5min <- rbind(A = colSums(flags_5min[A_sessions,, drop = FALSE]),
                    B = colSums(flags_5min[B_sessions,, drop = FALSE]))
unflagged_5min <- (208*2) - flags_5min
exclude_5min <- (unflagged_5min < 150) #exclude sessions with less than 150 volumes remaining

# Plot proportion of volumes scrubbed (T=10 min) -------------------------------
flags <- flags[c("test_LR1","test_LR2","test_RL1","test_RL2"),,,"first_5_mins",,1] 
sub_keep_inds <- !(dimnames(flags)[[2]] %in% as.character(sub_drop))
flags <- flags[,sub_keep_inds,,]

saveRDS(flags, file = file.path(dir_github,'results','5_flags','flags.rds'))

# Effective Sample Size --------------------------------------------------------

#read in estimated ACFs for each session
ACFs <- readRDS(file.path(dir_results,'results/4_FC_ACFs.rds'))

#reshape ACF arrays to match flag array
subj_ACFs <- gsub('_.+','',dimnames(ACFs[[1]])[[2]])
runs_ACFs <- paste0('test_',gsub('.+_','',dimnames(ACFs[[1]])[[2]]))
ACFs2 <- vector('list', length=3)
#loop over baseline denoising
for(bb in 1:nB){
  ACFs2[[bb]] <- array(c(ACFs[[bb]]),
                       dim=c(1185, 4, 1113),
                       dimnames = list(volume = 1:1185,
                                       run = unique(runs_ACFs),
                                       subject = unique(subj_ACFs)))
  subjects_keep <- as.character(aggFD$Subject[!aggFD$drop])
  ACFs2[[bb]] <- ACFs2[[bb]][,,subjects_keep]
  ACFs2[[bb]] <- aperm(ACFs2[[bb]], perm = c(2, 3, 1)) #permute dimensions to match flag array
}

#function to compute effective sample size from ACF
T_eff_fun <- function(acf, scrub){
  #construct TxT autocorrelation matrix
  cor_mat <- toeplitz(acf)
  #truncate and scrub
  nT <- max(which(!is.na(scrub)))
  cor_mat <- cor_mat[1:nT,1:nT]
  cor_mat <- cor_mat[!scrub[1:nT],!scrub[1:nT]]
  nT_scrub <- sum(!scrub[1:nT]) #number of vol remaining after scrubbing
  #compute effective sample size
  nT_scrub^2 / sum(cor_mat^2)
}

#compute effective sample size for each base method, scrubbing method, and run
dims <- dim(flags)[-3]
T_eff <- array(NA, dim = dims, dimnames = dimnames(flags)[-3])
names(dimnames(T_eff)) <- c('run','subject','scrub')
T_eff <- list(P32 = T_eff, P36 = T_eff, FIX = T_eff, FIX_GSR = T_eff)
for(d1 in 1:dims[1]){ #loop over runs
  for(d2 in 1:dims[2]){ #loop over subjects
    for(d3 in 1:dims[3]){ #loop over scrubbing levels
      scrub <- flags[d1,d2,,d3]
      for(bb in 1:nB){ #loop over base denoising methods
        acf <- ACFs2[[bb]][d1,d2,]
        T_eff[[bb]][d1, d2, d3] <- T_eff_fun(acf = acf, scrub = scrub)
      }
    }
  }
}

T_eff_df <- lapply(T_eff, as.data.frame.table, responseName = "T_eff")

#for baseline variance using nominal sample size -- number of time points contributing to the estimate
flags_nVol <- apply(!flags, c(1,2,4), sum, na.rm=TRUE) #sum over time points -- number of vol's retained
flags_nVol <- apply(flags_nVol, 2:3, sum, na.rm=TRUE) #sum over sessions
names(dimnames(flags_nVol)) <- c('subject', 'FD')
nVol_df <- as.data.frame.table(flags_nVol, responseName = 'nVol')
nVol_df$FD <- factor(nVol_df$FD, levels = FD_levels2)

baseNames <- baseNames[-1]

# Make data.frame of err values ------------------------------------------------
for (base in baseNames) {

  print(paste0('~~~~~~~~~~~~~~~~~~~~~~~~~~ ',base, ' ~~~~~~~~~~~~~~~~~~~~~~~~~~ '))

  #0. read in error values (est - ground truth), take abs value, and take weighted avg over test splits
  err <- readRDS(file.path(dir_err, paste0("withDVARS_err_", base, ".rds")))
  err_plus <- readRDS(file.path(dir_err, paste0("withDVARS_err_", base, "_plus", ".rds")))

  flags <- readRDS(file.path(dir_results, 'results/5_flags', 'withDVARS_flags.rds'))

  bias <- array(NA, dim(err$FC_err$A), dimnames = dimnames(err$FC_err$A))
  abs_err <- array(NA, dim(err$FC_err$A), dimnames = dimnames(err$FC_err$A))
  bias_plus <- array(NA, dim(err_plus$FC_err$A), dimnames = dimnames(err_plus$FC_err$A))
  abs_err_plus <- array(NA, dim(err_plus$FC_err$A), dimnames = dimnames(err_plus$FC_err$A))

  #compute weighted average absolute error and bias over A and B splits
  for (ss in seq(nScrub0)) { #loop over FD thresholds
    subsetName <- scrubNames0[ss]

    for (mm in seq(length(nT_seq_names))) {
      scanDur <- nT_seq_names[mm]
      unflagged <- !flags[,,,scanDur,subsetName,base,"FD"]
      unflagged_plus <- !flags[,,,scanDur,subsetName,base,"FDplus"]

      #weight given to A vs. B = total number of unflagged volumes
      A_weights <- apply(unflagged[A_sessions,,], 2, sum, na.rm=TRUE)
      B_weights <- apply(unflagged[B_sessions,,], 2, sum, na.rm=TRUE)
      A_weights_plus <- apply(unflagged_plus[A_sessions,,], 2, sum, na.rm=TRUE)
      B_weights_plus <- apply(unflagged_plus[B_sessions,,], 2, sum, na.rm=TRUE)

      #set excluded session weights to 0
      A_weights[exclude_5min["A",]] <- A_weights_plus[exclude_5min["A",]] <- 0
      B_weights[exclude_5min["B",]] <- B_weights_plus[exclude_5min["B",]] <- 0

      #compute denominator for weighted average
      total_weights <- A_weights + B_weights
      total_weights_plus <- A_weights_plus + B_weights_plus
      total_weights[total_weights == 0] <- NA
      total_weights_plus[total_weights_plus == 0] <- NA

      A_mat <- err$FC_err$A[,,scanDur,subsetName]
      B_mat <- err$FC_err$B[,,scanDur,subsetName]
      A_mat_plus <- err_plus$FC_err$A[,,scanDur,subsetName]
      B_mat_plus <- err_plus$FC_err$B[,,scanDur,subsetName]

      #compute weighted averages between A and B splits
      bias_mat <- ((A_mat * A_weights) + (B_mat * B_weights)) / total_weights
      abs_err_mat <- ((abs(A_mat) * A_weights) + (abs(B_mat) * B_weights)) / total_weights

      bias_mat_plus <- ((A_mat_plus * A_weights_plus) + (B_mat_plus * B_weights_plus)) / total_weights_plus
      abs_err_mat_plus <- ((abs(A_mat_plus) * A_weights_plus) + (abs(B_mat_plus) * B_weights_plus)) / total_weights_plus

      bias[,,scanDur,subsetName] <- bias_mat
      abs_err[,,scanDur,subsetName] <- abs_err_mat
      bias_plus[,,scanDur,subsetName] <- bias_mat_plus
      abs_err_plus[,,scanDur,subsetName] <- abs_err_mat_plus

      rm(A_mat, B_mat, A_mat_plus, B_mat_plus)
      rm(bias_mat, abs_err_mat, bias_mat_plus, abs_err_mat_plus)
    }
  }

  #1. combine err/err_plus and bias/bias_plus to include only FD2+, exclude DVARS only
  err_all <- abind::abind(abs_err, abs_err_plus[,,,5, drop=FALSE])
  bias_all <- abind::abind(bias, bias_plus[,,,5, drop=FALSE])
  dimnames(err_all)[[4]] <- dimnames(bias_all)[[4]] <- FD_levels2
  err_all2 <- err_all[,,,FD_levels2]
  err_all <- err_all[,,,FD_levels]
  bias_all <- bias_all[,,,FD_levels]
  rm(err, err_plus, abs_err, abs_err_plus, bias, bias_plus)

  #2. remove subjects in sub_drop
  sub_keep_inds <- !(dimnames(err_all)[[1]] %in% (as.character(sub_drop)))
  err_all <- err_all[sub_keep_inds,,,]
  err_all2 <- err_all2[sub_keep_inds,,,]
  bias_all <- bias_all[sub_keep_inds,,,]

  saveRDS(err_all, file.path(dir_results, 'results', paste0('err_all_', base, '.rds')))
  saveRDS(err_all2, file.path(dir_results, 'results', paste0('err_all2_', base, '.rds')))
  saveRDS(bias_all, file.path(dir_results, 'results', paste0('bias_all_', base, '.rds')))

  names(dimnames(err_all)) <- names(dimnames(err_all2)) <- names(dimnames(bias_all)) <- c('Subject','Edge','Duration','FD')
  dimnames(err_all)[[2]] <- dimnames(err_all2)[[2]] <- dimnames(bias_all)[[2]] <- 1:nEdge

  #3. summarize over edges/subjects

  # A) median/mean FC error (and bias) over edges (by subject)  -- for boxplots/tests and FC error vs duration line plots
  MAE_subj <- apply(err_all, c(1,3,4), median, na.rm = TRUE)
  rMSE_subj <- sqrt(apply(err_all^2, c(1,3,4), mean, na.rm = TRUE))
  bias_subj <- apply(bias_all, c(1,3,4), mean, na.rm = TRUE)

  # B) median and mean FC error over subjects (by edge) -- for matrix plots
  MAE_edge <- apply(err_all, 2:4, median, na.rm=TRUE) 
  inds_himo <- dimnames(err_all)[[1]] %in% sub_himo #indices of high-motion subjects
  MAE_edge_himo <- apply(err_all[inds_himo,,,], 2:4, median, na.rm=TRUE) 
  MAE_edge_lomo <- apply(err_all[!inds_himo,,,], 2:4, median, na.rm=TRUE) 
  rMSE_edge <- sqrt(apply(err_all^2, 2:4, mean, na.rm=TRUE)) 
  rMSE_edge_himo <- sqrt(apply((err_all[inds_himo,,,])^2, 2:4, mean, na.rm=TRUE)) 
  rMSE_edge_lomo <- sqrt(apply((err_all[!inds_himo,,,])^2, 2:4, mean, na.rm=TRUE)) 

  #4. convert arrays to data frames

  #median/mean over edges by subject
  MAE_subj <- as.data.frame.table(MAE_subj, responseName = 'MAE')
  rMSE_subj <- as.data.frame.table(rMSE_subj, responseName = 'rMSE')
  bias_subj <- as.data.frame.table(bias_subj, responseName = 'bias')
  MAE_subj$rMSE <- rMSE_subj$rMSE; rm(rMSE_subj)

  #median/mean over subjects by edges
  MAE_edge <- as.data.frame.table(MAE_edge, responseName = 'MAE'); MAE_edge$group <- 'all'
  MAE_edge_himo <- as.data.frame.table(MAE_edge_himo, responseName = 'MAE'); MAE_edge_himo$group <- 'himo'
  MAE_edge_lomo <- as.data.frame.table(MAE_edge_lomo, responseName = 'MAE'); MAE_edge_lomo$group <- 'lomo'
  MAE_edge <- rbind(MAE_edge, MAE_edge_himo, MAE_edge_lomo); rm(MAE_edge_himo, MAE_edge_lomo)
  rMSE_edge <- as.data.frame.table(rMSE_edge, responseName = 'rMSE'); rMSE_edge$group <- 'all'
  rMSE_edge_himo <- as.data.frame.table(rMSE_edge_himo, responseName = 'rMSE'); rMSE_edge_himo$group <- 'himo'
  rMSE_edge_lomo <- as.data.frame.table(rMSE_edge_lomo, responseName = 'rMSE'); rMSE_edge_lomo$group <- 'lomo'
  rMSE_edge <- rbind(rMSE_edge, rMSE_edge_himo, rMSE_edge_lomo); rm(rMSE_edge_himo, rMSE_edge_lomo)
  MAE_edge$rMSE <- rMSE_edge$rMSE; rm(rMSE_edge)

  #summarize over edges
  MAE_edge_summ <- MAE_edge %>% group_by(Duration, FD, group) %>%
    summarize(MAE = median(MAE, na.rm=TRUE),
              MSE = (mean(rMSE^2, na.rm=TRUE))) #for the mean does not matter if we summarize over edges or subjects first

  #save results for visualization
  save(MAE_edge, file = file.path(dir_results, 'results','7_MAE',paste0('MAE_edge_',base,'.RData')))
  save(MAE_subj, MAE_edge_summ, file = file.path(dir_github, 'results','7_MAE',paste0('MAE_',base,'.RData')))
  save(bias_subj, file = file.path(dir_github, 'results','7_MAE',paste0('bias_',base,'.RData')))

  #5. compute differences in error vs. no scrubbing 

  compute_diff <- function(df, FD_base, group_cols, value_col, pct=TRUE){
    if(is.character(value_col)) value_col <- which(names(df)==value_col)
    var <- names(df)[value_col]
    df <- df[,c(group_cols,value_col)] #remove extraneous variables
    df0 <- df[df$FD==FD_base,]
    var0 <- paste0(var,'0')
    names(df0)[names(df0)==var] <- var0
    df0$FD <- NULL
    df <- left_join(df, df0)
    vals <- unlist(df[,var])
    vals0 <- unlist(df[,var0])
    if(pct) return((vals - vals0)/vals0)
    if(!pct) return(vals - vals0)
  }

  MAE_subj$MAE_delta <- compute_diff(df = MAE_subj, FD_base = 'None', group_cols = 1:3, value_col = 'MAE')
  MAE_subj$rMSE_delta <- compute_diff(df = MAE_subj, FD_base = 'None', group_cols = 1:3, value_col = 'rMSE')
  MAE_edge$MAE_delta <- compute_diff(df = MAE_edge, FD_base = 'None', group_cols = c(1:3,5), value_col = 'MAE')
  MAE_edge$rMSE_delta <- compute_diff(df = MAE_edge, FD_base = 'None', group_cols = c(1:3,5), value_col = 'rMSE')

  #save results for visualization
  saveRDS(MAE_subj, file = file.path(dir_github, 'results','7_MAE',paste0('MAE_diff_subj_',base,'.rds')))
  saveRDS(MAE_edge, file = file.path(dir_results, 'results','7_MAE',paste0('MAE_diff_edge_',base,'.rds')))

  #6. compute BASELINE VARIANCE using nominal and effective sample size (T=10 min)

  # A) get ESS (nVol_eff) -- sum effective scan duration over runs, average over visits
  T_eff_sum <- T_eff_df[[base]] %>%
    group_by(subject, scrub) %>%
    summarize(nVol_eff = sum(T_eff)/2)
  names(T_eff_sum)[2] <- 'FD'
  save(T_eff_sum, file = file.path(dir_github, 'results','7_MAE',paste0('T_eff_sum_',base,'.RData')))

  # B) get squared error for each subject and edge
  errsq_10min <- as.data.frame.table((err_all2[,,"first_5_mins",])^2, responseName = 'sqerr'); rm(err_all2)

  # C) merge squared error with nVol and nVol_eff
  names(nVol_df)[1] <- names(T_eff_sum)[1] <- 'Subject'
  baselineSD <- left_join(errsq_10min, nVol_df); rm(errsq_10min) #merge in nominal sample size
  baselineSD <- left_join(baselineSD, T_eff_sum[,1:3]) #merge in effective sample size
  baselineSD$baselineSD <- sqrt(baselineSD$sqerr * baselineSD$nVol) # baseline variance = squared err * T
  baselineSD$baselineSD_eff <- sqrt(baselineSD$sqerr * baselineSD$nVol_eff) # baseline variance = squared err * T_eff

  # D) summarize over EDGES (mean/median) for each SUBJECT
  baselineSD$group <- ifelse(baselineSD$Subject %in% sub_himo, 'Above-Average Motion', 'Below-Average Motion')
  baselineSD$group <- factor(baselineSD$group, levels = c('Below-Average Motion','Above-Average Motion'))
  baselineSD_subj <- baselineSD %>%
    group_by(Subject, FD, group) %>%
    summarize(baselineSDnom_med = median(baselineSD, na.rm=TRUE),
              baselineSDeff_med = median(baselineSD_eff, na.rm=TRUE),
              baselineSDnom_mean = mean(baselineSD, na.rm=TRUE),
              baselineSDeff_mean = mean(baselineSD_eff, na.rm=TRUE))

  # E) summarize over SUBJECTS (mean/median) at each EDGE
  baselineSD_edge <- baselineSD %>%
    group_by(Edge, FD) %>%
    summarize(baselineSDnom_med = median(baselineSD, na.rm=TRUE),
              baselineSDeff_med = median(baselineSD_eff, na.rm=TRUE),
              baselineSDnom_mean = mean(baselineSD, na.rm=TRUE),
              baselineSDeff_mean = mean(baselineSD_eff, na.rm=TRUE))
  baselineSD_edge_grp <- baselineSD %>%
    group_by(Edge, FD, group) %>%
    summarize(baselineSDnom_med = median(baselineSD, na.rm=TRUE),
              baselineSDeff_med = median(baselineSD_eff, na.rm=TRUE),
              baselineSDnom_mean = mean(baselineSD, na.rm=TRUE),
              baselineSDeff_mean = mean(baselineSD_eff, na.rm=TRUE))
  baselineSD_edge$group <- 'all'
  baselineSD_edge <- baselineSD_edge[,c(1:2,7,3:6)]
  baselineSD_edge <- rbind(baselineSD_edge, baselineSD_edge_grp); rm(baselineSD_edge_grp)

  #7. compute differences in baseline SD vs. no scrubbing

  baselineSD_edge$baselineSDnom_med_delta <- compute_diff(df=baselineSD_edge, FD_base='None', group_cols=1:3, value_col=4)
  baselineSD_edge$baselineSDeff_med_delta <- compute_diff(df=baselineSD_edge, FD_base='None', group_cols=1:3, value_col=5)
  baselineSD_edge$baselineSDnom_mean_delta <- compute_diff(df=baselineSD_edge, FD_base='None', group_cols=1:3, value_col=6)
  baselineSD_edge$baselineSDeff_mean_delta <- compute_diff(df=baselineSD_edge, FD_base='None', group_cols=1:3, value_col=7)

  baselineSD_subj$baselineSDnom_med_delta <- compute_diff(df=baselineSD_subj, FD_base='None', group_cols=1:3, value_col=4)
  baselineSD_subj$baselineSDeff_med_delta <- compute_diff(df=baselineSD_subj, FD_base='None', group_cols=1:3, value_col=5)
  baselineSD_subj$baselineSDnom_mean_delta <- compute_diff(df=baselineSD_subj, FD_base='None', group_cols=1:3, value_col=6)
  baselineSD_subj$baselineSDeff_mean_delta <- compute_diff(df=baselineSD_subj, FD_base='None', group_cols=1:3, value_col=7)

  #save results for visualization
  save(baselineSD_subj, file = file.path(dir_github, 'results','7_MAE',paste0('baselineSD_',base,'.RData')))
  save(baselineSD_edge, file = file.path(dir_results, 'results','7_MAE',paste0('baselineSD_edge_',base,'.RData')))

  #8. determine duration to achieve "target" MSE at each edge/subject (based on lenient censoring at duration = 17.5 min)

  targetDur <- 17.5 #halfway between the min (5 min) and max (~30 min)
  targetDurText <- paste0('first_',targetDur/2,'_mins') #divide by 2 because total duration is split over 2 runs

  #grab target rMSE based on lenient censoring for each edge
  MAE_edge_target <- MAE_edge %>%
    filter(FD == 'FD5', Duration == targetDurText) %>%
    group_by(Edge, group) %>%
    select(Edge, group, rMSE_target = rMSE)
  MAE_edge <- MAE_edge %>%
    filter(FD != 'None') %>%
    left_join(MAE_edge_target) #bring in target MAE, defined at FD5 at T = 17.5 min

  #grab target rMSE based on lenient censoring for each subject
  MAE_subj <- MAE_subj %>% select(Subject, Duration, FD, rMSE)
  MAE_subj_target <- MAE_subj %>%
    filter(FD == 'FD5', Duration == targetDurText) %>%
    group_by(Subject) %>%
    select(Subject, rMSE_target = rMSE)
  MAE_subj <- MAE_subj %>%
    filter(FD != 'None') %>%
    left_join(MAE_subj_target) #bring in target MAE, defined at FD5 at T = 17.5 min

  # Function to identify the Duration that achieves target MAE

  #this function constructs a piecewise linear function between x and y
  #it then finds all the roots, i.e. the places where the function crosses y_target
  #it returns the minimum of all roots (the minimum duration required to achieve the target)
  #it will return 0 or Inf if the function does not cross the target
  find_min_root <- function(x, y, y_target) {
    if(all(y > y_target)) return(Inf)
    if(all(y < y_target)) return(0)
    roots <- numeric()
    for (i in 1:(length(x) - 1)) {
      x0 <- x[i]; x1 <- x[i + 1]
      y0 <- y[i]; y1 <- y[i + 1]
      # Check if the segment crosses y_target
      if ((y_target - y0) * (y_target - y1) <= 0 && y0 != y1) {
        # Linear interpolation to solve for x where f(x) = y_target
        t <- (y_target - y0) / (y1 - y0)
        x_root <- x0 + t * (x1 - x0)
        roots <- c(roots, x_root)
      }
    }
    return(min(roots))
  }

  format_duration <- function(x){ as.numeric(gsub("_mins", "", gsub("first_","",x)))*2 }
  MAE_edge$Duration <- format_duration(MAE_edge$Duration)
  MAE_subj$Duration <- format_duration(MAE_subj$Duration)

  #compute the minimum duration to achieve the target rMSE
  MAE_minDur_edge <- MAE_edge %>%
    group_by(Edge, FD, group) %>%
    summarize(target = rMSE_target[1], #all the same, just grab the first one
              minDur = find_min_root(Duration, rMSE, target))

  MAE_minDur_subj <- MAE_subj %>%
    group_by(Subject, FD) %>%
    summarize(target = rMSE_target[1], #all the same, just grab the first one
              minDur = find_min_root(Duration, rMSE, target))

  #what proportion of edges never achieve the target MAE?
  MAE_minDur_edge %>% group_by(FD) %>% summarize(mean(is.infinite(minDur)))

  #what proportion of edges require less than 5 minutes to achieve the target MAE?
  MAE_minDur_edge %>% group_by(FD) %>% summarize(mean((minDur==0)))

  ### Compute the percentage change in required duration (vs minDur for FD5)

  #first, deal with zeros and Infs
  MAE_minDur_edge$minDur[MAE_minDur_edge$minDur==0] <- 4 #if the minDur is too small to observe, set it to 4 min
  MAE_minDur_edge$minDur[is.infinite(MAE_minDur_edge$minDur)] <- 30 #if the minDur is too large to observe, set it to 30 min
  MAE_minDur_subj$minDur[MAE_minDur_subj$minDur==0] <- 4 #if the minDur is too small to observe, set it to 4 min
  MAE_minDur_subj$minDur[is.infinite(MAE_minDur_subj$minDur)] <- 30 #if the minDur is too large to observe, set it to 30 min

  #now, compute change in percent and minutes
  MAE_minDur_edge$Duration_change_pct <- compute_diff(df = MAE_minDur_edge, FD_base = 'FD5', group_cols = 1:3, value_col = 'minDur')
  MAE_minDur_edge$Duration_change_mins <- compute_diff(df = MAE_minDur_edge, FD_base = 'FD5', group_cols = 1:3, value_col = 'minDur', pct=FALSE)
  MAE_minDur_subj$Duration_change_pct <- compute_diff(df = MAE_minDur_subj, FD_base = 'FD5', group_cols = 1:2, value_col = 'minDur')
  MAE_minDur_subj$Duration_change_mins <- compute_diff(df = MAE_minDur_subj, FD_base = 'FD5', group_cols = 1:2, value_col = 'minDur', pct=FALSE)

  #what proportion of edges require longer duration vs. Lenient?
  MAE_minDur_edge %>% group_by(FD) %>% summarize(mean(Duration_change_mins > 0))

  #what proportion of edges require more than double the scan duration (> 100%)?
  MAE_minDur_edge %>% group_by(FD) %>% summarize(mean(Duration_change_pct > 1))

  #save for visualization
  save(MAE_minDur_subj, file = file.path(dir_github, 'results','7_MAE',paste0('minDuration_subj_',base,'.RData')))
  save(MAE_minDur_edge, file = file.path(dir_results, 'results','7_MAE',paste0('minDuration_edge_',base,'.RData')))

  #save a few edges for illustration of duration change
  set.seed(98765)
  edges_illustration <- sort(sample(1:nEdge, 3)) #18357 37994 61718
  MAE_edge_illustration <- MAE_edge %>% filter(Edge %in% edges_illustration, group == 'all')
  MAE_minDur_edge_illustration <- MAE_minDur_edge %>% filter(Edge %in% edges_illustration, group == 'all')
  save(edges_illustration, MAE_edge_illustration, MAE_minDur_edge_illustration, file = file.path(dir_github, 'results','7_MAE',paste0('minDurationIllustration_',base,'.RData')))

}

