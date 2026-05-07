####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Compute FC values
# How: First 1_MeanSignals.R, 1, ParcMatLabels.R, 2_ScrubMeas_DVARS.R, and 2_ScrubMeas_FD.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## -----------------------------------------------------------------------------
source("0_SharedCode.R")

for (with_S1200 in c(TRUE, FALSE)) {
  subjects <- if (with_S1200) { subjects_S1200 } else { subjects_RT }
  iters <- if (with_S1200) { iters_S1200 } else { iters_RT }
  visits <- if (with_S1200) { visits_S1200 } else { visits_RT }

  for (pproc in pprocs) {
    baseNames <- switch(pproc, MPP=MPP_baseNames, FIX=FIX_baseNames)
    nB <- length(baseNames)
    with_FIX <- pproc == "FIX"

    ## -------------------------------------------------------------------------

    flag_adjust_plus <- function(flag){
      stopifnot(is.logical(flag))
      nT <- length(flag)

      # Drop one back and two ahead.
      flag <- rbind(
        flag,
        c(flag[seq(2, nT)], FALSE),
        c(FALSE, flag[seq(nT-1)]),
        c(FALSE, FALSE, flag[seq(nT-2)])
      )
      flag <- colSums(flag)>0

      # Drop segments of less of than five volumes.
      fw <- which(flag)
      fw_dropidx <- which(diff(fw) < (5+1) & diff(fw) > 1)
      for (ff in fw_dropidx) {
        flag[seq(fw[ff], fw[ff+1])] <-  TRUE
      }

      flag
    }

    exp4 <- function(x) {
      x <- cbind(x, rbind(0, diff(x)))
      cbind(x, x^2)
    }


    ## -------------------------------------------------------------------------
    # Parcellation -------------------------------------------------------------
    ParcMat <- readRDS(parc_fname)
    cor_mask <- upper.tri(diag(parc_res2))

    ## -------------------------------------------------------------------------
    # Edit this function to control skipping iterations.
    skip_ii <- function(FC_fname) {
      # Will skip this iteration if the FC file(s) already exists.
      all(file.exists(FC_fname))
    }

    ## -------------------------------------------------------------------------
    get_FC <- function(this_flag, this_cii) {
      this_T <- length(this_flag)
      nScrub <- sum(this_flag)

      if (nScrub < this_T - 2) {
        return(cor(this_cii[!this_flag,])[cor_mask])
      } else {
        # Too many frames were scrubbed to calculate FC
        return(matrix(NA, nrow=parc_res2, ncol=parc_res2)[cor_mask])
      }
    }

    flag_wrap1_val <- function(my_flag, my_cii){
      list(FC=get_FC(my_flag, my_cii), flag=my_flag)
    }

    flag_wrap2_tmasks <- function(flag=NULL, nuis_bb=NULL, plus=FALSE, include_DVARS=FALSE) {
      if (is.null(flag)) {
        flag <- rep(FALSE, nT_1sess[length(nT_1sess)])
      } else {
        flag <- as.logical(flag)
        stopifnot(length(flag) == nT_1sess[length(nT_1sess)])
      }

      tmasks <- lapply(nT_seq, seq)

      flag <- lapply(tmasks, function(x){flag[x]})

      out <- vector("list", length=length(tmasks))
      names(out) <- names(tmasks)

      for (tt in seq(length(tmasks))) {
        vols_tt <- tmasks[[tt]]
        flag_tt <- flag[[tt]]
        T_tt <- length(flag_tt)

        if (include_DVARS) {
          stopifnot(all(names(tmasks) == names(dvars$dv)))
          dvflag_tt <- (dvars$dv[[tt]]$DPD > 5) & (dvars$dv[[tt]]$ZD > qnorm(1-.05/(nT_seq[tt])))
          flag_tt <- flag_tt | dvflag_tt
        }

        if (plus) { flag_tt <- flag_adjust_plus(flag_tt) }

        nScrub <- sum(flag_tt)

        #add in DCT regressors for temporal filtering
        nuis0_tt <- cbind(1, nT_DCT[[tt]])

        #for !FIX add in the other nuisance regressors
        if (!is.null(nuis_bb)) {
          nuis0_tt <- cbind(nuis0_tt, scale(nuis_bb[vols_tt,]))
        }

        # 107321_v2_LR: constant RPs. `scale` above makes `NA` values.
        nuis0_tt[,apply(is.na(nuis0_tt), 2, all)] <- 1

        #if (nScrub+ncol(nuis0_tt) < T_tt - 2) {
        if (nScrub < T_tt - 2) {
          if (nScrub > 0) {
            # One-hot encode outlier flags
            spikes <- matrix(0, nrow=T_tt, ncol=nScrub)
            spikes[seq(0, nScrub-1)*T_tt + which(flag_tt)] <- 1
          } else {
            # No scrubbing
            spikes <- NULL
          }

          nuis_tt <- cbind(nuis0_tt, spikes)
          cii2 <- nuisance_regression(cii[vols_tt,], nuis_tt)

        } else {
          # Too many volumes / DOF removed
          spikes <- NULL
          cii2 <- cii * NA
        }

        out[[tt]] <- flag_wrap1_val(flag_tt, cii2)
        #out[[tt]]$parc_ts <- cii2
      }
      out
    }

    ## -------------------------------------------------------------------------
    time <- Sys.time()

    fd_dir <- file.path(dir_scrubMeas, "FD")

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
      FC_fname <- file.path(dir_FC, baseNames, paste0(dvprefix, suffix, ".rds"))
      if (skip_ii(FC_fname)) { next }

      # CIFTI
      fname_prefix <- paste0("rfMRI_REST", visit, "_", acquisition)
      data_dir <- if (COMPUTER=="MPro") {
        file.path(subject, fname_prefix)
      } else {
        file.path(subject, "MNINonLinear", "Results", fname_prefix)
      }

      #Pre-computed regressors from 1_MeanSignals.R
      precomputed_fname <- paste0(subject, "_v", visit + (!test)*2, "_", acquisition, ".rds")
      ms <- try(readRDS(file.path(dir_meanSignals, precomputed_fname)))
      if (inherits(ms, "try-error")) { next }
      if (length(ms$cii) != 1200) { next } # 119732_v2_LR
      ms <- scale(as.data.frame(ms)[seq(nDrop+1, hcp_T),])

      if (with_FIX) {
        fnames <- list(
          CIFTI= file.path(data_dir, paste0(fname_prefix, "_Atlas_hp2000_clean.dtseries.nii"))
        )

        ms <- ms[,c("cii"),drop = FALSE] #global signal

      } else {
        fnames <- list(
          RP = file.path(data_dir, "Movement_Regressors.txt"),
          CIFTI = file.path(data_dir, paste0(fname_prefix, "_Atlas.dtseries.nii"))
        )

        ### Pre-computed regressors
        ms <- ms[,c("wm", "csf", "cii")] #mean white matter, mean csf, and global signal regressors
        cc <- try(readRDS(file.path(dir_CompCor, precomputed_fname))) #compcor regressors
        if (inherits(cc, "try-error")) { next }
        cc <- scale(cbind(cc$PCs$wm_cort[,seq(5)], cc$PCs$csf[,seq(5)], cc$PCs$wm_cblm[,seq(5)]))
      }

      # If retest, the data needs to be loaded.
      if (!test) {
        data_zip <- if (with_FIX) {
          file.path(
            dir_HCP_retest_archive,
            paste0(subject, "_3T_rfMRI_REST", visit, "_fixextended.zip")
          )
        } else {
          file.path(
            dir_HCP_retest_archive,
            paste0(subject, "_3T_rfMRI_REST", visit, "_preproc.zip")
          )
        }

        for (fname in fnames) {
          if (!file.exists(file.path(dir_HCP_retest, fname))) {
            cmd <- paste("unzip", data_zip, fname, "-d", dir_HCP_retest)
            system(cmd)
            stopifnot(file.exists(file.path(dir_HCP_retest, fname)))
          }
        }
      }

      #main directory containing the HCP data
      read_dir <- ifelse(test, dir_HCP_test, dir_HCP_retest)

      #read in nuisance parameters for non-FIX data
      if (!with_FIX) {
        rp <- try(as.matrix(read.table(file.path(read_dir, fnames$RP))))
        if (inherits(rp, "try-error")) { next }
        rp <- scale(cbind(rp, rp^2))
        rp <- scale(rp[seq(nDrop+1, hcp_T),])

        base_nuis <- list(
          P32 = cbind(rp, exp4(ms[,c("wm", "csf")])),
          P36 = cbind(rp, exp4(ms[,c("wm", "csf", "cii")]))
        )
        stopifnot(all(names(base_nuis) == baseNames))
      } else {
        base_nuis <- list(
          FIX = NULL,
          FIX_GSR = exp4(ms[,c("cii"), drop = FALSE])
        )
      }

      # Scrubbing measures
      cat("\tReading scrubbing measures.\n")
      fd <- try(readRDS(file.path(fd_dir, paste0("FD_", suffix, ".rds"))))
      if (inherits(fd, "try-error")) { next }
      if (with_DVARS) {
        dvars <- try(readRDS(file.path(
          dir_slate, "results/2_DVARS_MPP",
          paste0("DVARS_", subject, "_v", visit + (!test)*2,
                "_", acquisition, ".rds")
        )))
        if (inherits(dvars, "try-error")) { next }
      }

      # -------------------------------------------------------------------------------------------------------------------
      # Read CIFTI. Drop first 15 frames. Parcellate.
      # (Is equivalent to parcellating after nuisance regression.)
      if (!file.exists(file.path(dir_HCP_test, fnames$CIFTI))) { next }
      cat("\tReading CIFTI and parcellating.\n")
      cii <- try(
        as.matrix(read_xifti(
          file.path(dir_HCP_test, fnames$CIFTI), brainstructures="all"
        ))[,seq(nDrop+1, hcp_T)]
      )
      if (inherits(cii, "try-error")) { next }
      cii <- t(cii - rowMeans(cii)) %*% ParcMat

      for (bb in seq(nB)) {
        baseName_bb <- baseNames[bb]
        #nuis_bb <- if (with_FIX) { NULL } else { base_nuis[[bb]] }
        nuis_bb <- base_nuis[[bb]]
        FC_fname <- file.path(dir_FC, baseNames[bb], paste0(dvprefix, suffix, ".rds"))
        #if (file.exists(FC_fname)) { next }

        FCval_ii <- FCval_plus_ii <- vector("list")

        # Get FC values ----------------------------------------------------------
        # Nothing
        cat("\tBase")
        FCval_ii[["Base"]] <- FCval_plus_ii[["Base"]] <- flag_wrap2_tmasks(NULL, nuis_bb, plus=FALSE, include_DVARS=FALSE)

        if (with_DVARS) {
          FCval_ii[["DVARS only"]] <- flag_wrap2_tmasks(NULL, nuis_bb, plus=FALSE, include_DVARS=with_DVARS)
          FCval_plus_ii[["DVARS only"]] <- flag_wrap2_tmasks(NULL, nuis_bb, plus=TRUE, include_DVARS=with_DVARS)
        }

        # FD
        cat("\n\tFD")
        for (FD_cut in FD_cuts) {
          FD_type <- "og_nfc_l4"
          this_name <- paste("FD", FD_type, FD_cut, sep="___")
          flag <- fd[[FD_type]] > FD_cut
          FCval_ii[[this_name]] <- flag_wrap2_tmasks(flag, nuis_bb, plus=FALSE, include_DVARS=with_DVARS)
          FCval_plus_ii[[this_name]] <- flag_wrap2_tmasks(flag, nuis_bb, plus=TRUE, include_DVARS=with_DVARS)
        }

        saveRDS(FCval_ii, FC_fname)
        saveRDS(FCval_plus_ii, gsub(".rds", "_plus.rds", FC_fname, fixed=TRUE))
        print(length(unique(lapply(FCval_ii, '[[', "first_5_mins"))))
        cat("\n")
      }

      # Unload retest data.
      if (!test) {
        if (!with_S1200) { unlink(file.path(dir_HCP_retest, fnames$RP)) }
        unlink(file.path(dir_HCP_retest, fnames$CIFTI))
      }

      print(Sys.time() - time)
      time <- Sys.time()
      cat("\n")
    } # iters
  } # pprocs
} # with_S1200
