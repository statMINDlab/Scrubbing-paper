####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Compute DVARS
# How: No previous scripts needed
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## -----------------------------------------------------------------------------
source("0_SharedCode.R")

# loop over which dataset to analyze (with_S1200 = TRUE for full HCP, FALSE for 42 retest subjects)
for (with_S1200 in c(TRUE, FALSE)) {
  subjects <- if (with_S1200) { subjects_S1200 } else { subjects_RT }
  iters <- if (with_S1200) { iters_S1200 } else { iters_RT }
  visits <- if (with_S1200) { visits_S1200 } else { visits_RT }

  ## ---------------------------------------------------------------------------
  time <- Sys.time()

  for (ii in seq(nrow(iters))) {
    # Get iteration info -------------------------------------------------------
    subject <- iters[ii, "subject"]
    acquisition <- iters[ii, "acquisition"]
    test <- iters[ii, "test"]
    visit <- iters[ii, "visit"]

    # Get output file; skip if done --------------------------------------------
    dvars_fname <- file.path(
      dir_slate, "results/2_DVARS_MPP",
      paste0(
        "DVARS_", subject, "_v", visit + (!test)*2,
        "_", acquisition, ".rds"
      )
    )
    if (file.exists(dvars_fname)) { next }

    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    cat(paste0(
      "Subject ", subject, ", ", as.character(acquisition), " ",
      ifelse(test, "test", "retest"), " ", visit,
      " (", ii, " of ", nrow(iters), ")", "\n"
    ))
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    # Get input files. ---------------------------------------------------------
    # RPs
    fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
    data_dir <- switch(COMPUTER,
                       RED = file.path(subject, "MNINonLinear", "Results", fname_prefix),
                       MPro = file.path(subject, fname_prefix)
    )
    fnames <- list(
      CIFTI = file.path(data_dir, paste0(fname_prefix, "_Atlas.dtseries.nii"))
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

    if (with_S1200) { if (!file.exists(fnames$CIFTI)) { next } }

    # Read CIFTI. Drop first 15 frames. ------------------------------------------
    cat("\tReading data.\n")
    cii0 <- as.matrix(read_cifti(fnames$CIFTI, brainstructures="all"))
    if (ncol(cii0) != 1200) { next }
    cii0 <- t(cii0)[seq(nDrop+1, hcp_T),]
    cii_mean <- mean(cii0)

    # Compute DVARS for each session subset time. (The dual threshold will differ)
    get_DVARS <- function(nT_mm){
      cii <- t(cii0[seq(nT_mm),])
      cii <- cii - apply(cii, 1, mean) # for nuisance regression
      cii <- (cii + cii_mean) / cii_mean * 1000
      cii <- cii - apply(cii, 1, mean) # for DVARS normalization
      cii <- t(cii)
      DVARS(cii, normalize=FALSE)$measure
    }
    dv <- setNames(lapply(nT_seq, get_DVARS), nT_seq_names)
    saveRDS(list(dv=dv, grand=c(mean=cii_mean)), dvars_fname)

    # Unload retest data.
    if (!test) { unlink(fnames$CIFTI) }

    print(Sys.time() - time)
    time <- Sys.time()
  } # iters
} # with_S1200
