####################################################################
# Who: Amanda Mejia
# When: Fall 2024 to Spring 2025
# What: Obtain ground truth BWAS correlations & estimates with scrubbing and shorter duration
# How: First, 4_AggFCandFlags.R must be run to aggregate FC and flags
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

source("0_SharedCode.R")
setwd(dir_slate)

library(dplyr) #version 1.1.4

## ---------------------------------------------------------------------------

nS <- length(subjects_S1200)

## ---------------------------------------------------------------------------

# Part 0: Identify highly reliable behavioral measures

demo <- read.csv(file.path(dir_data_github,'unrestricted_HCP_demographics.csv'))
demo_RT <- read.csv(file.path(dir_data_github,'unrestricted_HCP_demographics_retest.csv'))
demo <- subset(demo, as.character(Subject) %in% as.character(demo_RT$Subject))
demo <- dplyr::arrange(demo, Subject)
demo_RT <- dplyr::arrange(demo_RT, Subject)

#redefine demographic/behavioral variables
y_vars_cat <- read.csv(file.path(dir_github, 'data', 'HCP_S1200_DataDictionary_Oct_30_2023.csv'))
y_vars_cat <- y_vars_cat[y_vars_cat$category %in% c("Cognition", "Emotion", "Motor", "Personality"), c("category", "columnHeader")]
y_vars_cat <- y_vars_cat[grep("Unadj|NEOFAC", y_vars_cat$columnHeader), ]
y_vars_cat <- as.character(y_vars_cat$columnHeader)

demo1 <- demo[,y_vars_cat]
demo2 <- demo_RT[,y_vars_cat]
var_tot <- (apply(demo1, 2, var, na.rm=TRUE) + apply(demo2, 2, var, na.rm=TRUE))/2
var_within <- apply(demo1 - demo2, 2, var, na.rm=TRUE) * (1/2)
ICC_demo <- 1 - var_within/var_tot
tail(sort(ICC_demo))
(yvar <- names(ICC_demo[which.max(ICC_demo)])) #CogTotalComp_Unadj
saveRDS(ICC_demo, file = file.path(dir_github, 'results', '8_BWAS', 'ICC_demo.rds'))

# Part 1: Get the ground truth BWAS rho value

# a) Get the behavioral measure for each subject: CogTotalComp_Unadj
# b) Get the ground truth FC for each subject: all 4 sessions (complete subjects only), FD2 scrubbing
# c) Compute cor(X, y). Account for age and sex.

# a) Get the behavioral measure for each subject: CogTotalComp_Unadj
demo_df <- read.csv(file.path(dir_data_github,'unrestricted_HCP_demographics.csv'))
demo_df <- demo_df[,c("Subject","Gender","Age",yvar)]
#reorder rows of demo_df to match subject ordering in FC array
row.names(demo_df) <- as.character(demo_df$Subject)
demo_df <- demo_df[subjects_S1200,]
y <- demo_df[,yvar]

for (bb in 1:nB) {

  baseName <- baseNames[bb]
  cat(baseName, "\n")

  df_rho <- NULL #initialize df to save rho vals
  gc()

  # b) Get the ground truth FC for each subject: all 4 sessions (complete subjects only), FD2 scrubbing

  #ground truth FC (use FD2)
  pfname_FD2 <- scrubNames[5] #partial file name
  fname_bs <- paste0('withS1200_withDVARS_FC_', baseName, '_', pfname_FD2, '.rds')
  FC_true_b <- readRDS(file.path(dir_results, 'results', '4_AggFC', fname_bs))[,,,1,1] #sess x sub x edge

  gc()
  X <- matrix(FC_true_b, nrow = 4)
  rm(FC_true_b); gc() #free up memory
  gc()

  # aside -- estimate ICC for T = 30 min
  X1 <- fishZ(matrix(colMeans(X[c(1,3),]), nrow = nS, ncol = nEdge)) #average over LR1 and RL1
  X2 <- fishZ(matrix(colMeans(X[c(2,4),]), nrow = nS, ncol = nEdge)) #average over LR2 and RL2
  var_tot <- (apply(X1, 2, var, na.rm=TRUE) + apply(X2, 2, var, na.rm=TRUE))/2
  var_within <- apply(X1 - X2, 2, var, na.rm=TRUE) * (1/2)
  ICC_X <- 1 - var_within/var_tot

  # project ICC estimate to T = 60 min
  var_sig <- var_tot - var_within #signal variance
  var_within_60min <- var_within/2 #double the scan duration --> half the noise variance
  var_tot_60min <- var_sig + var_within_60min #signal variance remains the same
  ICC_X_60min <- 1 - var_within_60min/var_tot_60min
  save(ICC_X, ICC_X_60min, file = file.path(dir_github, 'results', '8_BWAS', paste0('BWAS_',baseName,'_ICC_FC.RData')))

  #average over 4 sessions (subjects without 4 valid sessions will become NA)
  X <- matrix(colMeans(X), nrow = nS, ncol = nEdge)
  print(sum(is.na(X[,1])))
  rm(FC_true_b); gc() #free up memory
  X <- fishZ(X)

  # c) Compute cor(X, y). Do NOT adjust for age and sex so the correlations are as strong as possible so the GT is more precise

  #compute correlation between CogTotalComp_Unadj and FC
  keep <- (!is.na(y) & !is.na(rowSums(X))) #identify subjects to use in analysis

  print(sum(keep))

  rho_true <- apply(X, 2, function(x) cor(x[keep], y[keep]))

  ## ---------------------------------------------------------------------------
  # Part 2: Estimate rho with shorter scan duration

  durations <- c(2.5, 3.75, 5, 6.25, 7.5, 8.75, 10, 12.5, 14.22)
  nDur <- length(durations)

  #loop over scrubbing levels (None, FD5, FD2, FD2+)
  df_keep <- NULL
  for (ss in c(1,2,5,6)) {
    pfname_ss <- scrubNames[ss] #partial file name
    scrub_ss <- FD_levels2[ss]
    cat(scrub_ss, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    #loop over durations
    for(dd in 1:nDur){
      dur <- durations[dd]
      cat(dur, "minutes")
      pfname_sd <- gsub("14.22", dur, pfname_ss) #use first X min of LR and RL

      gc()
      fname_bs <- paste0('withS1200_withDVARS_FC_', baseName, '_', pfname_sd, '.rds')
      FC_bs <- readRDS(file.path(dir_results, 'results', '4_AggFC', fname_bs))[,,,1,1] #sess x sub x edge

      #repeat analysis for both visits (average error/bias values later)
      for(vv in 1:2){
        runs_v <- paste0(c("test_LR", "test_RL"), vv)
        #avg over runs (LR and RL) and Fisher transform
        X <- matrix(FC_bs[runs_v,,], nrow = 2);
        X <- matrix(colMeans(X), nrow = nS, ncol = nEdge)
        X <- fishZ(X)

        keep_sdv <- (!is.na(y) & !is.na(rowSums(X))) * keep
        print(sum(keep_sdv))
        rho_ss <- apply(X, 2, function(x) cor(x[keep], y[keep], use = 'complete'))
        df_rho_ss <- data.frame(base = baseName,
                                scrub = scrub_ss,
                                duration = dur,
                                visit = vv,
                                edge = 1:nEdge,
                                rho_est = rho_ss,
                                rho_true = rho_true)
        df_rho <- rbind(df_rho, df_rho_ss)
        saveRDS(df_rho, file = file.path(dir_slate, 'results', '8_BWAS', paste0('BWAS_',baseName,'.rds')))
        gc()

        df_keep_sdv <- data.frame(base = baseName,
                                  scrub = scrub_ss,
                                  duration = dur,
                                  visit = vv,
                                  nKeep = sum(keep_sdv))
        df_keep <- rbind(df_keep, df_keep_sdv)
        saveRDS(df_keep, file = file.path(dir_github, 'results', '8_BWAS', paste0('BWAS_',baseName,'_nKeep.rds')))
      }
      rm(FC_bs); gc() #free up memory

    } #end loop over durations
  } #end loop over scrubbing levels
} #end loop over baseline denoising methods

#split up large BWAS files to save to Github
for (bb in 1:nB) {

  baseName <- baseNames[bb]
  cat(baseName, "\n")

  df_rho <- readRDS(file.path(dir_slate, 'results', '8_BWAS', paste0('BWAS_',baseName,'.rds')))

  #loop over scrubbing levels, save subsets
  for (ss in c(1,2,5,6)) {

    pfname_ss <- scrubNames[ss] #partial file name
    scrub_ss <- FD_levels2[ss]
    cat(scrub_ss, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    df_rho_scrub <- subset(df_rho, scrub == scrub_ss)
    saveRDS(df_rho_scrub, file=file.path(dir_github, 'results', '8_BWAS', paste0('BWAS_',baseName,'_',scrub_ss,'.rds')))

  } #end loop over scrubbing levels
} #end loop over baseline denoising methods

#save ICC results for other scrub levels 
for (bb in 2:nB) {
  
  baseName <- baseNames[bb]
  cat(baseName, "\n")
  
  for (ss in c(1,2,5,6)) {
    pfname <- scrubNames[ss]
    scrub_ss <- FD_levels2[ss]
    cat(scrub_ss, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    fname_bs <- paste0('withS1200_withDVARS_FC_', baseName, '_', pfname, '.rds')
    FC_true_b <- readRDS(file.path(dir_results, 'results', '4_AggFC', fname_bs))[,,,1,1] #sess x sub x edge
    
    gc()
    X <- matrix(FC_true_b, nrow = 4)
    rm(FC_true_b); gc() #free up memory
    gc()
    
    # aside -- estimate ICC for T = 30 min
    X1 <- fishZ(matrix(colMeans(X[c(1,3),]), nrow = nS, ncol = nEdge)) #average over LR1 and RL1
    X2 <- fishZ(matrix(colMeans(X[c(2,4),]), nrow = nS, ncol = nEdge)) #average over LR2 and RL2
    var_tot <- (apply(X1, 2, var, na.rm=TRUE) + apply(X2, 2, var, na.rm=TRUE))/2
    var_within <- apply(X1 - X2, 2, var, na.rm=TRUE) * (1/2)
    ICC_X <- 1 - var_within/var_tot
    
    # project ICC estimate to T = 60 min
    var_sig <- var_tot - var_within #signal variance
    var_within_60min <- var_within/2 #double the scan duration --> half the noise variance
    var_tot_60min <- var_sig + var_within_60min #signal variance remains the same
    ICC_X_60min <- 1 - var_within_60min/var_tot_60min
    save(ICC_X, ICC_X_60min, file = file.path(dir_github, 'results', '8_BWAS', paste0('BWAS_',baseName,'_',scrub_ss,'_ICC_FC.RData')))
  }
}
