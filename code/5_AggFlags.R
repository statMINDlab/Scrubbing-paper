####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Compute aggregate flags to use in computing MAE and scrubbing rate
# How: First, 4_AggFCandFlags.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## ---------------------------------------------------------------------------
source("0_SharedCode.R")
setwd(dir_slate)

dvprefix <- "withDVARS_"

with_S1200 <- FALSE
subjects <- if (with_S1200) { subjects_S1200 } else { subjects_RT }
iters <- if (with_S1200) { iters_S1200 } else { iters_RT }
visits <- if (with_S1200) { visits_S1200 } else { visits_RT }
S1200_prefix <- if (with_S1200) { "withS1200_" } else { "" }

## ---------------------------------------------------------------------------
ParcMat <- readRDS(parc_fname)

## ---------------------------------------------------------------------------

aggFlags <- setNames(vector("list", nB), baseNames)
aggFlags <- list(FD = aggFlags, FDplus = aggFlags)

for (do_plus in c(FALSE, TRUE)) {
  plusName <- if (do_plus) { "_plus" } else { "" }

  for (bb in seq(nB)) {
    baseName <- baseNames[bb]

    #read in the aggregated flags for each baseline denoising method and scrubbing level
    aggFlags[[ifelse(!do_plus, "FD", "FDplus")]][[baseName]] <- setNames(lapply(scrubNames0, function(q){
      readRDS(file.path(dir_Agg, paste0(dvprefix, "flags_", baseName, "_", q, plusName, ".rds")))
    }), scrubNames0)
  }
}

#make a big array 
aggFlags <- lapply(aggFlags, function(q) { lapply(q, function(scrub) { abind::abind(scrub, along=5) }) }) #combine scrub levels
aggFlags <- lapply(aggFlags, function(base) { abind::abind(base, along=6) }) #combine base names
aggFlags <- abind::abind(aggFlags, along=7) 

#save the super-aggregated result
saveRDS(aggFlags, file.path(dir_slate, "results/5_flags", paste0(dvprefix, "flags.rds")))

## Plots comparing FD --------------------------------------------------------
library(gridExtra) #version 2.3

flags <- readRDS(file.path(dir_slate, 'results/5_flags/withDVARS_flags.rds'))

subjects <- dimnames(flags)[[2]]
sessions <- dimnames(flags)[[1]]

flag_summary <- data.frame(subject = character(),
                           session = character(),
                           n_expanded = numeric(),
                           n_stringent = numeric(),
                           ratio = numeric(),
                           stringsAsFactors = FALSE)

for (subject in subjects) {
  for (session in sessions) {
    f <- flags[session, subject, , "first_14.22_mins", "FD___og_nfc_l4___0.2", "P36", ]
    
    n_expanded <- sum(f[, "FDplus"])
    n_stringent <- sum(f[, "FD"])
    ratio <- ifelse(n_stringent == 0, NA, n_expanded / n_stringent)
    
    flag_summary <- rbind(flag_summary, data.frame(subject = subject,
                                                   session = session,
                                                   n_expanded = n_expanded,
                                                   n_stringent = n_stringent,
                                                   ratio = ratio))
  }
}

#remove rows with missing ratios
quantiles <- quantile(flag_summary$ratio, probs = c(0.05, 0.95), na.rm = TRUE)

#find sessions closest to each quantile
pick_sessions <- flag_summary[which.min(abs(flag_summary$ratio - quantiles[1])), ]
pick_sessions <- rbind(pick_sessions,
                       flag_summary[which.min(abs(flag_summary$ratio - quantiles[2])), ],
                       flag_summary[which.min(abs(flag_summary$ratio - quantiles[3])), ])

print(pick_sessions)

#create filtered_iters to get the correct ScrubMeas file 
sess <- pick_sessions$session
visit <- as.integer(sub(".*([12])$", "\\1", sess))       
test  <- grepl("^test_", sess)                           
acq   <- sub("^test_", "", sess)                        
acq   <- toupper(substr(acq, 1, 2))

filtered_iters <- data.frame(
  visit = visit,
  test = test,
  acquisition = factor(acq, levels = c("LR", "RL")), 
  subject = pick_sessions$subject,
  stringsAsFactors = FALSE
)
filtered_iters$suffix  <- paste0(filtered_iters$subject, "_v", filtered_iters$visit, "_", as.character(filtered_iters$acquisition))

FD_data <- list()

for (ii in seq(nrow(filtered_iters))){
  
  subject <- filtered_iters[ii, "subject"]
  acquisition <- filtered_iters[ii, "acquisition"]
  test <- filtered_iters[ii, "test"]
  visit <- filtered_iters[ii, "visit"]
  session <- paste0("test_", acquisition, visit)
  suffix <- filtered_iters[ii, "suffix"]

  fd_dir <- file.path(dir_scrubMeas, "FD")
  FD_ii <- readRDS(file.path(fd_dir, paste0("FD_", suffix, ".rds")))
  FD_og_ii <- FD_ii$og #original FD
  FD_lf_ii <- FD_ii$og_nfc_l4 #lagged and filtered FD
  
  flags_ii <- flags[session,subject,,"first_14.22_mins","FD___og_nfc_l4___0.2","P36",]
  FD2_ii <- flags_ii[,"FD"]
  FD2plus_ii <- flags_ii[,"FDplus"]
  
  FD_df_ii <- data.frame(FD_og = FD_og_ii,
                         FD_lf = FD_lf_ii,
                         FD2 = FD2_ii,
                         FD2plus = FD2plus_ii,
                         vol = 1:1185,
                         visit = visit,
                         acquisition = acquisition,
                         subject = subject)
  
  #new column indicating time points flagged with FD2+ (FD2plus = TRUE)
  FD_df_ii$FD2plus_vols <- ifelse(FD_df_ii$FD2plus == TRUE, FD_df_ii$vol, NA)
  
  #new column indicating time points where FD2 is also flagged 
  FD_df_ii$FD2_vols <- ifelse(FD_df_ii$FD2 == TRUE & FD_df_ii$FD2plus == TRUE, FD_df_ii$vol, NA)
  
  FD_data[[length(FD_data) + 1]] <- FD_df_ii
}

FD_df <- do.call(rbind, FD_data) 
FD_df <- FD_df[FD_df$FD_lf <= 1, ] 

FD_expanded_only <- FD_df[FD_df$FD2plus == TRUE & FD_df$FD2 == FALSE, ] #red dots, Expanded only 
FD_both <- FD_df[FD_df$FD2 == TRUE & FD_df$FD2plus == TRUE, ] #green dots, Expanded and Stringent

#plot comparing FD2 and FD2+
pdf(file.path(dir_github, 'plots', 'FD2_vs_FD2+.pdf'), width=8, height=5)
print(ggplot(FD_df, aes(x=vol,y=FD_lf)) +
  xlim(500,800) + 
  geom_hline(yintercept=0.2, col="darkgray") +
  geom_line(linewidth=0.4, color = 'grey40') + 
  geom_point(data=FD_expanded_only, aes(x=vol, y=FD_lf, fill="#b63545"), size=1.25, shape=4, color="#b63545", stroke=0.8) + 
  geom_point(data=FD_both, aes(x=vol, y=FD_lf, fill="#b3ce91"), size=1.25, shape=21, color="#b3ce91", stroke=0.3) + 
  scale_fill_manual(values = c("#b63545" = "#b63545", "#b3ce91" = "#b3ce91")) +
  coord_cartesian(ylim = c(0, 0.7)) +
  theme_few() + 
  facet_wrap(~subject, ncol=1) + 
  labs(x='Volume', y='Filtered FD', title='FD2 and FD2+') + 
  theme(legend.position='none'))
dev.off()

#plot comparing original and lagged & filtered FD
pdf(file.path(dir_github, "plots", "FDog_vs_FDlf.pdf"), width=8, height=5)
print(ggplot(FD_df, aes(x=vol)) +
        xlim(400,600) +
        geom_line(aes(y=FD_og), color='gray') +
        geom_line(aes(y=FD_lf), color='gray20') +
        coord_cartesian(ylim = c(0, 0.7)) +
        theme_few() +
        facet_wrap(~subject, ncol=1) +
        labs(x="Volume", y="FD", title="Original vs. Lagged/Filtered FD"))
dev.off()
