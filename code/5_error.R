####################################################################
# Who: Damon Pham, Joanne Hwang
# When: Spring 2024 - Fall 2025
# What: Compute error measures
# How: First, 4_AggFCandFlags.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## -----------------------------------------------------------------------------
source("0_SharedCode.R")
setwd(dir_slate)

## -----------------------------------------------------------------------------
ParcMat <- readRDS(parc_fname)

## -----------------------------------------------------------------------------

n_FC <- (parc_res2^2 - parc_res2)/2

time <- Sys.time()

weights <- TRUE
baseNames <- baseNames[-1] # P36, FIX, FIX_GSR

# Handle subjects

for (do_plus in c(FALSE, TRUE)) {
  plusName <- if (do_plus) { "_plus" } else { "" }

  for (bb in seq(nB)) {
    baseName <- baseNames[bb]
    cat(baseName, "\n\n")

    err_fname <- file.path(dir_err, paste0(dvprefix, "err_", baseName, plusName, ".rds"))
    #if (file.exists(err_fname)) { next }

    # Get the retest estimates to use for the ground truth
    # four retest sessions x  42 subjects x 90K connections
    FC_gt <- readRDS(file.path(dir_Agg, paste0(dvprefix, "FC_", baseName, "_FD___og_nfc_l4___0.2_first14.22mins.rds")))
    FC_gt[] <- psych::fisherz(FC_gt)
    stopifnot(all(dimnames(FC_gt)[[1]] == visits_RT))
    visits_gt <- list(
      A = grepl("retest", visits_RT) | (visits_RT %in% c("test_LR2", "test_RL2")), #retest + visit 2 of main HCP
      B = grepl("retest", visits_RT) | (visits_RT %in% c("test_LR1", "test_RL1")) #retest + visit 1 of main HCP
    )
    FC_gt <- list(A = FC_gt[visits_gt$A,,,,], #retest + visit 2 of main HCP
                  B = FC_gt[visits_gt$B,,,,]) #retest + visit 1 of main HCP
    stopifnot(all(vapply(FC_gt, function(q){dim(q)[1]==6}, FALSE)))

    # average over sessions, ignoring any bad sessions
    FC_gt_unweighted <- lapply(FC_gt, colMeans, na.rm=TRUE)

    if (weights) {
      # get weights (# unflagged volumes per session)
      flags <- readRDS(file.path(dir_slate, 'results/5_flags', 'withDVARS_flags.rds'))
      unflagged <- !flags[,,,"first_14.22_mins","FD___og_nfc_l4___0.2",baseName,"FD"]
      A_sessions <- dimnames(FC_gt$A)[[1]] #which runs
      B_sessions <- dimnames(FC_gt$B)[[1]] #which runs
      subjects <- dimnames(unflagged)[[2]]

      unflagged_A <- unflagged[A_sessions,,]
      unflagged_B <- unflagged[B_sessions,,]

      A_weights <- matrix(NA, length(A_sessions),length(subjects),dimnames=list(A_sessions,subjects))
      B_weights <- matrix(NA, length(B_sessions),length(subjects),dimnames=list(B_sessions,subjects))
      A_weights <- apply(unflagged_A, c(1,2), sum, na.rm=FALSE) #keep session and subject names
      B_weights <- apply(unflagged_B, c(1,2), sum, na.rm=FALSE)

      # average over weighted sessions
      FC_gt_weighted <- list(A = matrix(NA, dim(FC_gt$A)[2], dim(FC_gt$A)[3], dimnames=list(dimnames(A_weights)[[2]])),
                             B = matrix(NA, dim(FC_gt$B)[2], dim(FC_gt$B)[3], dimnames=list(dimnames(B_weights)[[2]])))

      for (subject in 1:dim(FC_gt$A)[2]) { #num subjects is same for A and B
        subject_data_A <- FC_gt$A[,subject,]
        subject_data_B <- FC_gt$B[,subject,]
        subject_weights_A <- A_weights[,subject]
        subject_weights_B <- B_weights[,subject]
        FC_gt_weighted$A[subject,] <- apply(subject_data_A, 2, weighted.mean, subject_weights_A, na.rm=TRUE)
        FC_gt_weighted$B[subject,] <- apply(subject_data_B, 2, weighted.mean, subject_weights_B, na.rm=TRUE)
      }
    }
    #end of computing ground truth

    FC_gt <- if (weights) { FC_gt_weighted } else { FC_gt_unweighted }

    if (!do_plus) {
      # average ground truth FC for matrices
      FC_combined <- abind::abind(FC_gt$A, FC_gt$B, along = 3) #combine A and B
      avg_FC_combined <- apply(FC_combined, c(1,2), mean) #average over A and B 
      FC_over_subj <- apply(avg_FC_combined, 2, mean) #average over subjects for visualization 
      
      abs_FC_over_subj <- abs(FC_over_subj)
      print(mean(abs_FC_over_subj))
      print(median(abs_FC_over_subj))
      
      saveRDS(FC_over_subj, file.path(dir_github, 'results', '5_error', paste0('FCgt_over_subj_',baseName,'.rds')))
    }
    
    flags <- readRDS(file.path(dir_slate, 'results/5_flags', 'withDVARS_flags.rds'))
    A_sessions <- c("test_LR1","test_RL1")
    B_sessions <- c("test_LR2","test_RL2")

    FC_err <- setNames(vector("list", nScrub0), scrubNames0)
    for (ss in seq(nScrub0)) {
      subsetName <- scrubNames0[ss]
      cat(subsetName, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

      agg_fname <- file.path(dir_Agg, paste0(dvprefix, "FC_", baseName, "_", subsetName, "_first2.5mins", plusName, ".rds"))
      if (!file.exists(agg_fname)) { next }

      FC_err[[ss]] <- setNames(vector("list", length(nT_seq_names)), nT_seq_names)
      for (mm in seq(length(nT_seq_names))) {
        cat("\t", nT_seq_names[mm], "\n")
        FC <- readRDS(gsub("first2.5mins", gsub("_", "", nT_seq_names[mm]), agg_fname))
        FC[] <- psych::fisherz(FC[])

        # get weights (# unflagged volumes per session)
        FD <- if (do_plus) { "FDplus" } else { "FD" }
        unflagged <- !flags[,,,nT_seq_names[mm],subsetName,baseName,FD]

        unflagged_A <- unflagged[A_sessions,,]
        unflagged_B <- unflagged[B_sessions,,]
        A_weights <- apply(unflagged_A, c(1,2), sum, na.rm=TRUE) 
        B_weights <- apply(unflagged_B, c(1,2), sum, na.rm=TRUE)

        FC_weighted <- list(A = matrix(NA, dim(FC)[2], dim(FC)[3], dimnames = list(dimnames(FC)[[2]])),
                            B = matrix(NA, dim(FC)[2], dim(FC)[3], dimnames = list(dimnames(FC)[[2]])))

        for (subject in 1:dim(FC)[2]) { #num subjects is same for A and B
          subject_data_A <- FC[A_sessions,subject,,,]
          subject_data_B <- FC[B_sessions,subject,,,]
          subject_weights_A <- A_weights[,subject]
          subject_weights_B <- B_weights[,subject]
          FC_weighted$A[subject,] <- apply(subject_data_A, 2, weighted.mean, subject_weights_A, na.rm=TRUE)
          FC_weighted$B[subject,] <- apply(subject_data_B, 2, weighted.mean, subject_weights_B, na.rm=TRUE)

          rm(subject_data_A, subject_data_B, subject_weights_A, subject_weights_B)
        }

        rm(unflagged, unflagged_A, unflagged_B, A_weights, B_weights)

        FC$A <- (FC_weighted$A - FC_gt$A)
        FC$B <- (FC_weighted$B - FC_gt$B)
        FC_err[[ss]][[mm]] <- FC
      }

      FC_err[[ss]] <- list(
        A = abind::abind(lapply(FC_err[[ss]], '[[', "A"), along=3),
        B = abind::abind(lapply(FC_err[[ss]], '[[', "B"), along=3)
      )
    }

    FC_err <- list(
      A = abind::abind(lapply(FC_err, '[[', "A"), along=4),
      B = abind::abind(lapply(FC_err, '[[', "B"), along=4)
    )
    saveRDS(list(FC_err = FC_err, FC_gt = FC_gt), err_fname)
  }
}
