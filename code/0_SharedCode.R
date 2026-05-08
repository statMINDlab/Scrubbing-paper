####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Defines variables/directories needed in later scripts
# How: Run this first before running other scripts
# Where: Can be run on HPC or locally, since it only requires access to Github files
####################################################################

# This is not necessary if you open the R project (*.Rproj) in RStudio
setwd('~/Documents/Github/Scrubbing-paper')

# Change these depending on where files are in your computer -------------------
COMPUTER <- ifelse(
  grepl("Quartz", getwd()), "RED", #HPC system
  ifelse(grepl("Documents", getwd()), "Personal") #your local computer
)

secret_fnames <- "0_FilePaths.txt" #edit this according to your system
if (!file.exists(secret_fnames)) {
  if (file.exists(file.path("code", secret_fnames))) {
    secret_fnames <- file.path("code", secret_fnames)
  } else if (file.exists(file.path("..", secret_fnames))) {
    secret_fnames <- file.path("..", secret_fnames)
  }
}
secret_fnames <- readLines(secret_fnames)

# HCP data, for HPC
dir_HCP_test <- switch(COMPUTER,
                       RED = secret_fnames[1], #HPC system
                       Personal = NULL)

dir_HCP_retest_archive <- secret_fnames[2]
# where to unzip HCP retest files for temporary use
dir_HCP_retest <- secret_fnames[3]

# Connectome Workbench
wb_path <- switch(COMPUTER,
                  RED = secret_fnames[4],
                  Personal = secret_fnames[5])


dir_github <- '~/Documents/Github/Scrubbing-paper/'
dir_slate <- '/N/project/LevScrubbing/Scrubbing-Cost'
dir_results <- dir_slate
dir_results_Github <- file.path(dir_github, 'results')
dir_project_ScrubPham2023 <- "/N/project/LevScrubbing" #location of results from an earlier project (Pham et al. 2023, Less is More) to be reused here


# Below should be the same -----------------------------------------------------
dir_data_github <- file.path(dir_github, "data") #in Github for local access
#dir_results_misc <- file.path(dir_github, "results") #in Github for local access
dir_Parc <- file.path(dir_data_github, "Parc") #in Github
dir_CompCor <- file.path(dir_project_ScrubPham2023, "analysis-results/1_CompCor") #on cluster storage
dir_meanSignals <- file.path(dir_project_ScrubPham2023, "analysis-results/1_MeanSignals") #on cluster storage
dir_scrubMeas <- file.path(dir_project_ScrubPham2023, "analysis-results/2_ScrubMeas") #on cluster storage
dir_FC <- file.path(dir_slate, "results/3_FC") # on cluster storage
dir_Agg <- file.path(dir_slate, "results/4_AggFC") # on cluster storage
dir_err <- file.path(dir_slate, "results/5_error") # on cluster storage
dir_plots <- file.path(dir_github, "plots") #on Github

parc_cii_fname <- file.path(dir_Parc, "Schaefer2018_400Parcels_Kong2022_17Networks_order.dlabel.nii")
hcp_dg_fname <- file.path(dir_data_github, "data", "data/unrestricted_HCP_demographics.csv") # gitignore

# The 45 Retest set subjects
# subjects_RT <- list.files(dir_HCP_test)
subjects_RT <- c(
  103818, 187547,
  105923, 192439,
  111312, 194140,
  114823, 195041,
  115320, 200109,
  122317, 200614,
  125525, 204521,
  130518, 250427,
  135528, 287248,
  137128, 341834,
  139839, 433839,
  143325, 562345,
  144226, 599671,
  146129, 601127,
  149337, 627549,
  149741, 660951,
  151526, 662551,
  158035, 783462,
  169343, 859671,
  172332, 861456,
  175439, 877168,
  177746, 917255,
  185442
)
omit_subjects_RT <- c(
  627549, # does not have retest zip files
  341834, # does not have retest FIX files in the zip (REST1 RL)
  143325  # retest MPP CIFTI for REST2_RL is truncated to 939 timepoints
)
subjects_RT <- subjects_RT[!(subjects_RT %in% omit_subjects_RT)]
rm(omit_subjects_RT)

# The ~1200 non-retest subjects
subjects_S1200_fname <- file.path(dir_data_github, "HCP_S1200_subjects.rds")
if (!file.exists(subjects_S1200_fname)) {
  subjects_S1200 <- sort(list.files(dir_HCP_test))
  subjects_S1200 <- subjects_S1200[!is.na(as.numeric(subjects_S1200))]
  saveRDS(subjects_S1200, subjects_S1200_fname)
}
subjects_S1200 <- readRDS(subjects_S1200_fname)

# Packages
library(fMRIscrub) #version 0.14.5
stopifnot(utils::packageVersion("fMRIscrub") >= "0.14.5")
library(ciftiTools) #version 0.16.1
stopifnot(utils::packageVersion("ciftiTools") >= "0.14.0")
ciftiTools.setOption("wb_path", wb_path)
library(fMRItools) #version 0.4.7
library(RColorBrewer) #version 1.1.3
library(ggplot2) #version 4.0.1
library(scales) #version 1.4.0
library(ggthemes) #version 5.1.0

hcp_TR <- .72
hcp_T <- 1200
nDrop <- 15

# `iters`
iters_RT <- expand.grid(
  visit=seq(2),
  test=c(TRUE, FALSE),
  acquisition=c("LR", "RL"),
  subject=subjects_RT
)
iters_S1200 <- expand.grid(
  visit=seq(2),
  acquisition=c("LR", "RL"),
  subject=subjects_S1200,
  test=TRUE
)

# Names of visits.
iters0 <- expand.grid( # temp
  visit=seq(2), test=c('test', 'retest'), acquisition=c("LR", "RL")
)
visits_RT <- paste0(iters0$test,'_',iters0$acquisition, iters0$visit)
iters0_S1200 <- expand.grid( # temp
  visit=seq(2), test='test', acquisition=c("LR", "RL")
)
visits_S1200 <- paste0(iters0_S1200$test,'_',iters0_S1200$acquisition, iters0_S1200$visit)
rm(iters0, iters0_S1200)

# Schaefer parcellation (400 parcels)
parc_fname <- file.path(dir_Parc, "ParcMat.rds")
parc_res <- 400
parc_res2 <- parc_res + 19 # plus 19 subcortical regions
nEdge <- (parc_res2^2 - parc_res2)/2  #87571

pprocs <- c("MPP", "FIX")
MPP_baseNames <- c("P32", "P36")
MPP_baseNames_long <- c("32P","36P")
FIX_baseNames <- c("FIX", "FIX_GSR") #add in GSR
FIX_baseNames_long <- c("ICA-FIX","FIX + GSR")
baseNames <- c(MPP_baseNames, FIX_baseNames)
baseNames_long <- c(MPP_baseNames_long, FIX_baseNames_long)
nB <- length(baseNames)

#always use DVARS along with FD for scrubbing
with_DVARS <- TRUE
dvprefix <- "withDVARS_"

FD_cuts <- rev(c(.2, .3, .4, .5))
FD_cuts_int <- as.integer(FD_cuts * 10)
FD_levels <- c('None','FD5','FD2','FD2+')
FD_levels_long <- c('None','Lenient', 'Stringent', 'Expanded')
FD_levels2 <- c('None','FD5','FD4','FD3','FD2','FD2+')
FD_levels2_long <- c('None','FD5\n(lenient)', 'FD4', 'FD3', 'FD2\n(stringent)', 'FD2+\n(stringent+\nexpanded)')

#file paths (within downstream scripts can replace "first14.22mins")
scrubNames0 <- c("Base", paste0("FD___og_nfc_l4___0.", FD_cuts_int)) #without FD2+ and "first14.22mins"
scrubNames <- paste0(scrubNames0, "_first14.22mins")
scrubNames <- c(scrubNames,
                paste0(scrubNames[length(scrubNames)], "_plus")) #add FD2+
#scrubNames2 <- gsub("___og_nfc_l4___0.", "", scrubNames, fixed=TRUE)
nScrub <- length(scrubNames)
nScrub0 <- length(scrubNames0) #5 (Base, FD5, FD4, FD3, FD2)

# Fisher transformation
fishZ <- function(r){ 0.5*log((1+r)/(1-r)) }

# Define session subset times
nT_1sess <- hcp_T - nDrop # 1200 - 15 == 1185 (drop first 15 frames)
nT_seq_mins <- c(2.5, 3.75, 5, 6.25, 7.5, 8.75, 10, 12.5, 14.22)
nT_seq_names <- paste0("first_", nT_seq_mins, "_mins")
nT_seq <- setNames(round(nT_seq_mins*60/hcp_TR), nT_seq_names)
stopifnot(nT_seq[length(nT_seq)] == nT_1sess)
nT_DCT <- lapply(nT_seq, function(q){
  fMRItools::dct_bases(q, fMRItools::dct_convert(q, hcp_TR, f=.01))
})

myCols <- c("black", "#5585bd", "#b3ce91", "#b63545") #for None / FD5 / FD2 / FD2+
myCols2 <- c("black", "deepskyblue", "blue", "purple", "orange","red")#for None / FD5 / FD4 / FD3 / FD2 / FD2+

