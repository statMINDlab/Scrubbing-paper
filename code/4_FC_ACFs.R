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

  #where the FC will be written

  for (bb in seq(nB)) {

    ts_fname <- file.path(dir_FC, baseNames[bb], paste0(dvprefix, suffix, "_ts.rds"))
    if (!file.exists(file.path(dir_FC, baseNames[bb], paste0(dvprefix, suffix, "_ts.rds")))) { next }
    parc_ts_ii <- readRDS(file=ts_fname) #T x (num_edge) matrix

    #estimating the autocorrelation function
    # acfs_ii <- apply(parc_ts_ii, 2, acf, plot=FALSE, lag.max=1000)
    # acfs_ii <- sapply(acfs_ii, function(x) c(x$acf)) #1001 x 419
    # acfs_ii <- apply(acfs_ii, 1, mean)
    acfs_ii <- apply(parc_ts_ii, 2, acf, plot=FALSE, lag.max=100) #truncate lags after 100
    acfs_ii <- sapply(acfs_ii, function(x) c(x$acf)) #101 x 419
    acfs_ii <- apply(acfs_ii, 1, median)

    ACFs[[bb]][1:101,ii] <- acfs_ii

    # # AR values
    # ar_ii <- apply(parc_ts_ii, 2, ar, order.max=10, aic=FALSE)
    # ar_ii <- sapply(ar_ii, function(x) c(x$ar)) #10 x 419
    # ar_ii <- apply(ar_ii, 1, mean)
    #
    # ARs[[bb]][1:10,ii] <- ar_ii

    # # Robust ACF (nearly identical to non-robust)
    # rob_acfs_ii <- apply(parc_ts_ii, 2, acfrob, lag.max=100, plot=FALSE)
    # rob_acfs_ii <- sapply(rob_acfs_ii, function(x) c(x$acf))
    # rob_acfs_ii <- apply(rob_acfs_ii, 1, median)
    #
    # rob_ACFs[[bb]][1:101,ii] <- rob_acfs_ii

  } #end loop over baseNames
} #end loop over subject

saveRDS(ACFs, file = file.path(dir_slate, "results", "4_FC_ACFs.rds"))
