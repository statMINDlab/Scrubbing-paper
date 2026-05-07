####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Compute FD scrubbing measures
# How: No previous scripts needed
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## ---------------------------------------------------------------------------
source("0_SharedCode.R")

# Check for correct version of gsignal
stopifnot(packageVersion("gsignal") <= "0.3-6.9000")

# Download an older version of gsignal, if needed
#install.packages("devtools")
#devtools::install_github("gjmvanboxtel/gsignal@937a29758d8c11c5f5f52a3c9cd047fa261a7d49") # 20220619 v0.3.6

# loop over which dataset to analyze (with_S1200 = TRUE for full HCP, FALSE for 42 retest subjects)
for (with_S1200 in c(TRUE, FALSE)) {
  subjects <- if (with_S1200) { subjects_S1200 } else { subjects_RT }
  iters <- if (with_S1200) { iters_S1200 } else { iters_RT }
  visits <- if (with_S1200) { visits_S1200 } else { visits_RT }

## ---------------------------------------------------------------------------

  nfc <- c(.31, .43) * .72 * 2
  nfiltc <- gsignal::cheby2(2, Rs=20, w=nfc, type="stop")

  time <- Sys.time()

  for (ii in seq(nrow(iters))) {
    # Get iteration info -------------------------------------------
    subject <- iters[ii, "subject"]
    acquisition <- iters[ii, "acquisition"]
    test <- iters[ii, "test"]
    visit <- iters[ii, "visit"]

    # Get output file; skip if done --------------------------------
    fd_fname <- file.path(
      dir_scrubMeas, "FD",
      paste0(
        "FD_", subject, "_v", visit + (!test)*2,
        "_", acquisition, ".rds"
      )
    )
    if (file.exists(fd_fname)) { next }

    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    cat(paste0(
      "Subject ", subject, ", ", as.character(acquisition), " ",
      ifelse(test, "test", "retest"), " ", visit,
      " (", ii, " of ", nrow(iters), ")", "\n"
    ))
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    # Get input files. ---------------------------------------------
    # RPs
    fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
    data_dir <- switch(COMPUTER,
                       RED = file.path(subject, "MNINonLinear", "Results", fname_prefix),
                       MPro = file.path(subject, fname_prefix)
    )
    fnames <- list(
      RP = file.path(data_dir, "Movement_Regressors.txt"),
      RP_dt = file.path(data_dir, "Movement_Regressors_dt.txt")
    )

    # If retest, the data needs to be loaded.
    if (!test) {
      data_zip_MPP <- file.path(
        dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip")
      )

      for (fname in fnames) {
        if (!file.exists(file.path(dir_HCP_retest, fname))) {
          cmd <- paste("unzip", data_zip_MPP, fname, "-d", dir_HCP_retest)
          system(cmd)
          stopifnot(file.exists(file.path(dir_HCP_retest, fname)))
        }
      }
    }

    fnames <- lapply(fnames, function(x){file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), x)})
    if (with_S1200) { if (!file.exists(fnames$RP)) { cat(fnames$RP, ": No RPs\n"); next } }

    rp <- as.matrix(read.table(fnames$RP)[,seq(6)])
    rp_nfc <- gsignal::filtfilt(nfiltc, rp)

    fd <- list(
      og = fMRIscrub::FD(rp),
      og_nfc_l4 = fMRIscrub::FD(rp_nfc, lag=4)
    )

    fd <- lapply(fd, function(x){x$measure[seq(nDrop+1, hcp_T)]})

    # Save FD.
    saveRDS(fd, fd_fname)

    # Unload retest data.
    if (!test) {
      unlink(fnames$RP)
      unlink(fnames$RP_dt)
    }

    print(Sys.time() - time)
    time <- Sys.time()
  }
}
