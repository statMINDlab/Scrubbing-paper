####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Create large FC array and FD array
# How: First, 3_FC.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## ---------------------------------------------------------------------------
source("0_SharedCode.R")
setwd(dir_slate)

writeFlags <- TRUE
writeFC <- TRUE

dvprefix <- "withDVARS_"

for (with_S1200 in c(TRUE, FALSE)) {
  subjects <- if (with_S1200) { subjects_S1200 } else { subjects_RT }
  iters <- if (with_S1200) { iters_S1200 } else { iters_RT }
  visits <- if (with_S1200) { visits_S1200 } else { visits_RT }
  
  S1200_prefix <- if (with_S1200) { "withS1200_" } else { "" }

  ## ---------------------------------------------------------------------------
  ParcMat <- readRDS(parc_fname)

  ## ---------------------------------------------------------------------------
  nEdge <- (parc_res2^2 - parc_res2)/2

  for (do_plus in c(FALSE, TRUE)) {
    plusName <- if (do_plus) { "_plus" } else { "" }

    #loop over baseline denoising methods
    for (bb in seq(nB)) {
      baseName <- baseNames[bb]
      cat(baseName, "\n\n")
      
      nScrub0 <- if (with_S1200 & do_plus) { 1 } else { 5 }
      
      #loop over scrubbing levels
      for (ss in seq(nScrub0)) {
        scrubName <- if (with_S1200 & do_plus) { scrubNames0[5] } else { scrubNames0[ss] }
        cat(scrubName, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        
        for (ii in seq(nrow(iters))) {
          # Get iters info
          subject <- iters[ii, "subject"]
          acquisition <- as.character(iters[ii, "acquisition"])
          test <- iters[ii, "test"]
          visit <- iters[ii, "visit"]
          
          # Read in FC values & scrubbing flags
          FC_fname <- paste0(dvprefix, subject, "_v", visit + (!test)*2, "_", acquisition, plusName, ".rds")
          FC_fname <- file.path(dir_FC, baseName, FC_fname)
          cat("\t", basename(FC_fname), "\n")
          #cat("\t", visits_RT[ifelse(ii%%8 == 0, 8, ii %%8)], "\n")
          if (!file.exists(FC_fname)) { warning("Missing file: ", FC_fname); next }
          FC <- readRDS(FC_fname) # a nested list: for each base denoising method and duration, contains FC values and flag
          
          # Initialize big FC array, if first iteration
          if (ii == 1) {
            methods <- names(FC)[grepl(scrubName, names(FC))]
            n_methods <- length(methods)
            print("n_methods:"); print(n_methods)
            if(writeFC){ FCagg <- array(NA, c(nrow(iters), nEdge, length(nT_seq), n_methods)) }
            if(writeFlags) { flags <- array(NA, c(nrow(iters), nT_1sess, length(nT_seq), n_methods)) }
          }
          
          if(writeFC){
            q <- lapply(FC[methods], function(x){do.call(cbind, lapply(x, '[[', "FC"))})
            q <- abind::abind(q, along=3)
            FCagg[ii,,,] <- q
          }
          
          if (writeFlags) {
            q <- lapply(FC[methods], function(x){do.call(cbind, lapply(
              x, function(y){c(y$flag, rep(NA, nT_1sess - length(y$flag)))}))})
            q <- abind::abind(q, along=3)
            flags[ii,,,] <- q
          }
        }
        
        cat("Saving...\n")
        
        if(writeFC) {
          dim(FCagg) <- c(length(visits), length(subjects), nEdge, length(nT_seq), n_methods)
          dimnames(FCagg) <- list(visits, subjects, NULL, nT_seq_names, methods)
          # Split FC by session subset length, for smaller file sizes
          FCagg_fname <- file.path(dir_Agg, paste0(S1200_prefix, dvprefix, "FC_", baseName, "_", scrubName, "_XX", plusName, ".rds"))
          for (mm in seq(length(nT_seq_names))) {
            mm_name <- gsub("_", "", nT_seq_names[mm])
            saveRDS(FCagg[,,,mm,,drop=FALSE], gsub("XX", mm_name, FCagg_fname))
          }
        }
        
        if (writeFlags) {
          dim(flags) <- c(length(visits), length(subjects), nT_1sess, length(nT_seq), n_methods)
          dimnames(flags) <- list(visits, subjects, NULL, nT_seq_names, methods)
          flags_fname <- file.path(dir_Agg, paste0(S1200_prefix, dvprefix, "flags_", baseName, "_", scrubName, plusName, ".rds"))
          saveRDS(flags, flags_fname)
        }
      } #end loop over scrubbing levels
    } #end loop over baseline denoising methods
  } #end loop over expanded scrubbing ("plus")
}
