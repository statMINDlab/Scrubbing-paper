####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Compute mean signals data
# How: No previous scripts needed
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## -----------------------------------------------------------------------------
source("0_SharedCode.R")

for (with_S1200 in c(TRUE, FALSE)) {
  subjects <- if (with_S1200) { subjects_S1200 } else { subjects_RT }
  iters <- if (with_S1200) { iters_S1200 } else { iters_RT }
  visits <- if (with_S1200) { visits_S1200 } else { visits_RT }

## -----------------------------------------------------------------------------
  time <- Sys.time()

  for (ii in seq(nrow(iters))) {
    subject <- iters[ii, "subject"]
    acquisition <- as.character(iters[ii, "acquisition"])
    test <- iters[ii, "test"]
    visit <- iters[ii, "visit"]

    # Get output file.
    ms_fname <- file.path(
      dir_meanSignals,
      paste0(subject, "_v", visit + (!test)*2, "_", acquisition, ".rds")
    )
    if (file.exists(ms_fname)) { next }

    cat(paste0(
      "Subject ", subject, ", ", as.character(acquisition), " ",
      ifelse(test, "test", "retest"), " ", visit,
      " (", ii, " of ", nrow(iters), ")", "\n"
    ))

    # Get input files.
    fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
    data_dir <- file.path(subject, "MNINonLinear", "Results", fname_prefix)
    fnames <- list(
      CIFTI = file.path(data_dir, paste0(fname_prefix, "_Atlas.dtseries.nii")),
      CIFTI_FIX = file.path(data_dir, paste0(fname_prefix, "_Atlas_hp2000_clean.dtseries.nii")),
      NIFTI = file.path(data_dir, paste0(fname_prefix, ".nii.gz")),
      NIFTI_labs = file.path(data_dir, "../../ROIs/Atlas_wmparc.2.nii.gz")
    )

    # If retest, the data needs to be loaded.
    if (!test) {
      
      # data_zip_MPP <- file.path(
      #   dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip")
      # )

      for (fname in fnames[c("CIFTI", "CIFTI_FIX", "NIFTI")]) {
        
        #CIFTI_FIX is in a different place
        data_zip_MPP <- if (fname == fnames$CIFTI_FIX) { file.path(dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_fixextended.zip"))
          } else {
          file.path(dir_HCP_retest_archive, paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip"))
          }
        
        if (!file.exists(file.path(dir_HCP_retest, fname))) {
          cmd <- paste("unzip", data_zip_MPP, fname, "-d", dir_HCP_retest)
          system(cmd)
          stopifnot(file.exists(file.path(dir_HCP_retest, fname)))
        }
      }
    }

    nii <- file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$NIFTI)
    nii <- try(RNifti::readNifti(nii))
    if (inherits(nii, "try-error")) { next }
    niiLabs <- file.path(dir_HCP_test, fnames$NIFTI_labs)
    niiLabs <- try(RNifti::readNifti(niiLabs))
    if (inherits(niiLabs, "try-error")) { next }
    noiseROIs <- fMRItools:::get_NIFTI_ROI_masks(file.path(dir_HCP_test, fnames$NIFTI_labs), c("wm_cort", "csf", "wm_cblm"))
    noiseROIs <- list(wm=noiseROIs$wm_cort | noiseROIs$wm_cblm, csf=noiseROIs$csf, wholebrain=niiLabs>0)
    noise_erosion <- list(wm=2, csf=1, wholebrain=0)
    for (rr in seq(length(noiseROIs))) {
      noiseROIs[[rr]][,,] <- as.logical(noiseROIs[[rr]]) * 1
      noiseROIs[[rr]][,,] <- fMRItools::erode_mask_vol(noiseROIs[[rr]], noise_erosion[[rr]], c(-1, 0, NA, NaN))
      noiseROIs[[rr]] <- apply(matrix(nii[noiseROIs[[rr]] > 0], ncol=1200), 2, mean)
    }
    rm(nii)

    #compute global signal
    cii <- file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$CIFTI)
    cii <- tryCatch({read_xifti(cii, brainstructures="all", flat=TRUE)}, error=function(cond){NA})
    if (identical(cii, NA)) { cat("\tBad CIFTI"); stop }
    noiseROIs$cii <- colMeans(as.matrix(cii))

    #compute global signal for FIX data
    cii_FIX <- file.path(ifelse(test, dir_HCP_test, dir_HCP_retest), fnames$CIFTI_FIX)
    cii_FIX <- tryCatch({read_xifti(cii_FIX, brainstructures="all", flat=TRUE)}, error=function(cond){NA})
    if (identical(cii_FIX, NA)) { cat("\tBad CIFTI"); stop }
    noiseROIs$cii_FIX <- colMeans(as.matrix(cii_FIX))

    saveRDS(noiseROIs, ms_fname)
    rm(cii)
    rm(cii_FIX)

    # Unload retest data.
    if (!test) {
      unlink(file.path(dir_HCP_retest, fnames$NIFTI))
      unlink(file.path(dir_HCP_retest, fnames$CIFTI))
      unlink(file.path(dir_HCP_retest, fnames$CIFTI_FIX))
    }

    print(Sys.time() - time)
    time <- Sys.time()
  }
}
