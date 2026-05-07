####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Compute aggregate FD data
# How: First, 2_ScrubMeas_FD.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## -----------------------------------------------------------------------------
source("0_SharedCode.R")

for (with_S1200 in c(TRUE, FALSE)) {
  subjects <- if (with_S1200) { subjects_S1200 } else { subjects_RT }
  iters <- if (with_S1200) { iters_S1200 } else { iters_RT }
  visits <- if (with_S1200) { visits_S1200 } else { visits_RT }

  S1200_prefix <- if (with_S1200) { "withS1200_" } else { "" }

  fd_dir <- file.path(dir_scrubMeas, "FD")
  fd_agg <- vector("list", length(iters))

  for (ii in seq(nrow(iters))) {
    # Get iteration info -------------------------------------------------------
    subject <- iters[ii, "subject"]
    acquisition <- as.character(iters[ii, "acquisition"])
    test <- iters[ii, "test"]
    visit <- iters[ii, "visit"]
    suffix <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition)
    cat(suffix, "\n")
    fd_agg[[ii]] <- try(readRDS(
      file.path(fd_dir, paste0("FD_", suffix, ".rds"))
    ))
    if (inherits(fd_agg[[ii]], "try-error")) {
      fd_agg[[ii]] <- list(og_nfc_l4=rep(NA, nT_1sess))
    }
  }
  fd_agg <- lapply(fd_agg, '[[', "og_nfc_l4")
  fd_agg <- do.call(rbind, fd_agg)
  saveRDS(fd_agg, file.path(dir_slate, "results/3_AggFD",
                            paste0(S1200_prefix, "AggFD.rds")))
}
