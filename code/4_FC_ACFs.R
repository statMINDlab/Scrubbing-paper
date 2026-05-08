####################################################################
# Who: Joanne Hwang and Amanda Mejia
# When: Fall 2024
# What: Compute ACFs to use in computing MAE
# How: First, 3_FC_timeseries.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## -----------------------------------------------------------------------------
source("0_SharedCode.R")

with_S1200 <- TRUE
subjects <- if (with_S1200) { subjects_S1200 } else { subjects_RT }
iters <- if (with_S1200) { iters_S1200 } else { iters_RT }
visits <- if (with_S1200) { visits_S1200 } else { visits_RT }

# initialize ACFs
ACFs <- matrix(0, nrow=1185, ncol=nrow(iters))
colnames(ACFs) <- paste0(iters$subject,'_',iters$acquisition,iters$visit)
ACFs <- rep(list(ACFs), nB)
names(ACFs) <- baseNames

for (ii in seq(nrow(iters))) {
  # Get iteration info -------------------------------------------
  subject <- iters[ii, "subject"]
  acquisition <- as.character(iters[ii, "acquisition"])
  test <- iters[ii, "test"]
  visit <- iters[ii, "visit"]

  cat(paste0(
    "Subject ", subject, ", ", acquisition, " ",
    ifelse(test, "test", "retest"), " ", visit,
    " (", ii, " of ", nrow(iters), ")", "\n"
  ))

  # Get files; skip if done --------------------------------
  suffix <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition)

  for (bb in seq(nB)) {

    ts_fname <- file.path(dir_FC, baseNames[bb], paste0(dvprefix, suffix, "_ts.rds"))
    if (!file.exists(file.path(dir_FC, baseNames[bb], paste0(dvprefix, suffix, "_ts.rds")))) { next }
    parc_ts_ii <- readRDS(file=ts_fname)

    #estimating the autocorrelation function
    acfs_ii <- apply(parc_ts_ii, 2, acf, plot=FALSE, lag.max=100) #truncate lags after 100
    acfs_ii <- sapply(acfs_ii, function(x) c(x$acf))
    acfs_ii <- apply(acfs_ii, 1, median)

    ACFs[[bb]][1:101,ii] <- acfs_ii
    
  } 
} 

saveRDS(ACFs, file = file.path(dir_slate, "results", "4_FC_ACFs.rds"))
