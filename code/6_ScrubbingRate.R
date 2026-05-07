####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Determine average scrubbing rate and compare data after scrubbing for each scrubbing method
# How: First, 3_AggFD.R and 5_AggFlags.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## ---------------------------------------------------------------------------
source("0_SharedCode.R")
setwd(dir_slate)

library(dplyr)

with_S1200 <- FALSE
subjects <- if (with_S1200) { subjects_S1200 } else { subjects_RT }
iters <- if (with_S1200) { iters_S1200 } else { iters_RT }
visits <- if (with_S1200) { visits_S1200 } else { visits_RT }
S1200_prefix <- if (with_S1200) { "withS1200_" } else { "" }

# Read in data. -----
### `aggFD` -----
# session by volume matrix of (mod)FD values.
aggFD <- readRDS(file.path(dir_results, "results/3_AggFD/AggFD.rds"))
# add summary stats (across volume) of `aggFD` to `iters`.
iters$FD_mean <- rowMeans(aggFD)
iters$FD2_rate <- rowMeans(aggFD > .2) #for each session, percent of volumes scrubbed
# `FDstats`: take the mean of the summary stats across sessions, for each subject.
FDstats <- aggregate(iters[c("FD_mean", "FD2_rate")], list(iters$subject), mean)
colnames(FDstats)[1] <- "Subject"

### `flags` -----
# `flags:` T/F
# visits * subjects * time * time subsets * scrubbing * base nuisance regression * plus
flags <- readRDS(file.path(dir_results, "results/5_flags" , paste0(dvprefix, "flags.rds")))
names(dimnames(flags)) <- c("visit", "subject", "volume", "time_subset", "FD_cutoff", "Baseline", "FD_type")

# Check that `flags` (from FC results) vs. `aggFD` (from before FC calc) is as expected.
q1 <- aggFD > .2
q2 <- matrix(flags[,,,9,"FD___og_nfc_l4___0.2","P32","FD"], nrow=336)
# q2 includes DVARS, q1 does not.
stopifnot(all(q2 - q1 >= 0))

# `nScrub`: number of scrubbed in this session, on average. list by time subset. then sum.
nScrub <- lapply(
  setNames(seq(dim(flags)[4]), nT_seq_names),
  function(q){apply(flags[,,seq(nT_seq[q]),q,,,], c(1,2,4,5,6), sum)}
)
nScrub <- abind::abind(nScrub, along=6)
stopifnot(!any(is.na(nScrub)))
nScrub <- aperm(nScrub, c(1,2,6,3,4,5))

##### Checks
# Check and update: 0 scrubbed for base (no scrubbing). remove/reorder this dim.
if (!with_DVARS) {
  stopifnot(all(nScrub[,,,"Base",,] == 0))
  nScrub <- nScrub[,,,dimnames(nScrub)[[4]]!="Base",,]
} else {
  nM <- dim(nScrub)[4]
  stopifnot(all(dimnames(nScrub)[[4]][seq(1)] == c("Base")))
  nScrub <- nScrub[,,,c(seq(2,nM), 1),,]
}

# Check: FDplus always flags more than FD
stopifnot(all(nScrub[,,,,,"FDplus"] - nScrub[,,,,,"FD"]) >= 0)

# Check and update: baseline has no effect on number scrubbed
stopifnot(length(unique(
  lapply(c("P32", "P36","FIX"), function(z){nScrub[,,,,z,]})
)) == 1)
# So we can drop that dimension
nScrub <- nScrub[,,,,1,]

# Check: increasing scrubbing threshold always flags more
stopifnot(all(nScrub[,,,"FD___og_nfc_l4___0.2",] - nScrub[,,,"FD___og_nfc_l4___0.3",]) >= 0)
stopifnot(all(nScrub[,,,"FD___og_nfc_l4___0.3",] - nScrub[,,,"FD___og_nfc_l4___0.4",]) >= 0)
stopifnot(all(nScrub[,,,"FD___og_nfc_l4___0.4",] - nScrub[,,,"FD___og_nfc_l4___0.5",]) >= 0)

# Ground truth & dropping subjects ----------

# Identify test LR 1 and test RL 1 sessions
idx_LR1 <- 1 + seq(0, 42-1)*8
idx_RL1 <- 5 + seq(0, 42-1)*8
# Checks
stopifnot(nrow(unique(iters[idx_LR1,c("visit", "test", "acquisition")])) == 1)
stopifnot(all(iters[idx_LR1,"subject"] == subjects))
stopifnot(nrow(unique(iters[idx_RL1,c("visit", "test", "acquisition")])) == 1)
stopifnot(all(iters[idx_RL1,"subject"] == subjects))

### Two ground-truth partitions of sessions. -----
visits_gt <- list(
  A = grepl("retest", visits_RT) | (visits_RT %in% c("test_LR2", "test_RL2")),
  B = grepl("retest", visits_RT) | (visits_RT %in% c("test_LR1", "test_RL1"))
)

# Array of number scrubbed
q <- lapply(visits_gt, function (z){
  colSums(nScrub[z,,"first_14.22_mins","FD___og_nfc_l4___0.2","FD"]) }) #using FD2 ground truth
# Time remaining per subject after scrubbing
nT_postScrub_gt <- do.call(rbind, lapply(q, function(z){ ((nT_1sess*6)-z)*hcp_TR/60 }))

### Drop subjects with insufficient data -----
# (require >50 minutes in both GT partitionings)
sub_drop <- nT_postScrub_gt[,apply(nT_postScrub_gt < 50, 2, any), drop = FALSE]
colnames(sub_drop) # 103818, 151526, 175439

# Drop subjects with insufficient data from `FDstats`
FDstats$drop <- FDstats$Subject %in% colnames(sub_drop)

### High vs low movers. -----
# Two options: meanFD and ScrubRate.  Based on Motion.pdf plot, ScrubRate seems to distinguish better.
med_meanFD <- median(subset(FDstats, !drop)$FD_mean) #0.1367949
avg_meanFD <- mean(subset(FDstats, !drop)$FD_mean) #0.1455963
#FDstats$motion <- (FDstats$FD_mean > med_meanFD)
# Defined by rate of scrubbing with FD2.
med_ScrubRate <- median(subset(FDstats, !drop)$FD2_rate) #0.1593354
avg_ScrubRate <- mean(subset(FDstats, !drop)$FD2_rate) #0.179242
FDstats$motion <- (FDstats$FD2_rate > avg_ScrubRate)

FDstats$motion <- factor(
  ifelse(FDstats$drop, "Drop", ifelse(FDstats$motion, "High Motion", "Low Motion")),
  levels=c("Low Motion", "High Motion", "Drop")
)
# FDstats$motion <- factor(FDstats$motion,
#                          levels = c("Low Motion", "High Motion", "Drop"),
#                          labels = c('Below-Average Motion', 'Above-Average Motion', 'Drop'))
# FDstats$motion2 <- factor(
#   ifelse(FDstats$drop, "Drop", ifelse(FDstats$motion2, "High Motion", "Low Motion")),
#   levels=c("Low Motion", "High Motion", "Drop")
# )
print(table(FDstats$motion))
#  Low Motion High Motion        Drop
#          21          18           3

FDstats <- arrange(FDstats, FD2_rate)
FDstats$Subject <- factor(FDstats$Subject, levels=FDstats$Subject)
pdf(file.path(dir_plots, "MotionSplit.pdf"), height=4, width=5)
ggplot(subset(FDstats, !drop), aes(x=Subject, y = FD2_rate)) +
  geom_bar(stat='identity', aes(fill=motion)) +
  ylab('% Volumes with FD > 0.2mm') + xlab('Subject ID') +
  scale_fill_manual('', values = c('gray','gold'), labels = c('Below-Average Motion', 'Above-Average Motion')) +
  theme_few() + theme(axis.text.x = element_text(angle=90, size=8), legend.position = 'bottom')
ggplot(subset(FDstats, !drop), aes(x=Subject, y = FD_mean)) +
  geom_bar(stat='identity', aes(fill=motion)) +
  ylab('Mean FD') + xlab('Subject ID') +
  scale_fill_manual('', values = c('gray','gold'), labels = c('Below-Average Motion', 'Above-Average Motion')) +
  theme_few() + theme(axis.text.x = element_text(angle=90, size=8), legend.position = 'bottom')
dev.off()

FDstats <- FDstats[match(subjects, FDstats$Subject),] #
saveRDS(FDstats, file.path(dir_slate, "results/6_scrubRate/AggFD_withSubjectSplit.rds"))
saveRDS(FDstats, file.path(dir_github, "results/6_scrubRate/AggFD_withSubjectSplit.rds"))

#mean FD scrubbing rate, and exclusion flag for each HCP retest subject
sub_drop <- as.numeric(as.character(FDstats$Subject[FDstats$drop])) #these 4 subjects have too little data for "ground truth"
sub_himo <- as.numeric(as.character(FDstats$Subject[FDstats$motion == "High Motion"]))
sub_lomo <- as.numeric(as.character(FDstats$Subject[FDstats$motion == "Low Motion"]))
save(sub_drop, sub_himo, sub_lomo, file = file.path(dir_github,'data','subjects_retest.RData'))

# Plots

### Time remaining in GT partitions per subject -----

# Make `arr_gtr`, an array of time remaining per subject in the 2 ground truth
#   partitions.
q <- lapply(visits_gt, function(z){ colSums(nScrub[z,as.character(FDstats$Subject),,,]) })
q$A <- abind::abind(A=q$A, along=5)
q$B <- abind::abind(B=q$B, along=5)
q <- abind::abind(q$A, q$B, along=5)
# sub_alt_order <- names(sort(rowMeans(
#   q[,"first_14.22_mins","FD___og_nfc_l4___0.2","FDplus",])))
# time (min.) remaining after scrubbing
arr_gtr <- (q*(-1) + rep(nT_seq*6, each=prod(dim(q)[seq(1)])) )*hcp_TR/60
arr_gtr <- arr_gtr[,"first_14.22_mins",,,]
names(dimnames(arr_gtr)) <- c("Subject", "FD_cutoff", "FD_type", "Sess_split")

# Convert `arr_gtr` to a data.frame.
q_df <- as.data.frame.table(arr_gtr, responseName = "Time")
levels(q_df$FD_cutoff) <- gsub(".*_", "", levels(q_df$FD_cutoff))

sub_by_meanFD <- FDstats$Subject[order(FDstats$FD_mean)]
sub_by_FDrate <- FDstats$Subject[order(FDstats$FD2_rate)]

color_m <- c(FD5 = "deepskyblue", 
             FD4 = "blue", 
             FD3 = "purple", 
             FD2 = "orange",
             Base = "black")

color_gg <- setNames(color_m[c("FD5","FD4","FD3","FD2","Base")], NULL)

pdf(file.path(dir_plots, "ScrubbingRate.pdf"), width=9, height=5)

q_df1 <- q_df
q_df1$Subject <- factor(q_df1$Subject, levels=sub_by_meanFD)
plt1 <- ggplot(q_df1, aes(x=Subject, y=Time)) +
  geom_point(aes(color=FD_cutoff, shape=FD_type)) +
  geom_hline(yintercept=30, lty=2) +
  geom_vline(xintercept = 1:42+0.5, color='lightgray') +
  scale_color_manual(values = color_gg) +
  scale_shape_manual(values = c(16, 4)) +
  ylab('Minutes remaining after scrubbing') + xlab('Subjects') +
  scale_y_continuous(limits=c(0,nT_1sess*6*hcp_TR/60), expand=c(0,0)) +
  ggthemes::theme_few() + facet_grid(.~Sess_split) +
  theme(axis.text.x = element_text(angle = 45, size=6, hjust=1)) +
  ggtitle("Time remaining for ground truth estimate after scrubbing\n(Subjects sorted by mean FD)")

q_df2 <- q_df
q_df2$Subject <- factor(q_df2$Subject, levels=sub_by_FDrate) # set subject order
plt2 <- ggplot(q_df2, aes(x=Subject, y=Time)) +
  geom_point(aes(color=FD_cutoff, shape=FD_type)) +
  geom_hline(yintercept=30, lty=2) +
  geom_vline(xintercept = 1:42+0.5, color='lightgray') +
  scale_color_manual(values = color_gg) +
  scale_shape_manual(values = c(16, 4)) +
  ylab('Minutes remaining after scrubbing') + xlab('Subjects') +
  scale_y_continuous(limits=c(0,nT_1sess*6*hcp_TR/60), expand=c(0,0)) +
  ggthemes::theme_few() + facet_grid(.~Sess_split) +
  theme(axis.text.x = element_text(angle = 45, size=6, hjust=1)) +
  ggtitle("Time remaining for ground truth estimate after scrubbing\n(Subjects sorted by FD flag rate)")

print(plt1)
print(plt2)

dev.off()

pdf(file.path(dir_plots, "Motion.pdf"), width=10, height=5)

q_df1 <- subset(FDstats, !drop)
q_df1$Subject <- factor(q_df1$Subject, levels=rev(sub_by_meanFD))
plt1 <- ggplot(q_df1, aes(x=Subject, y=FD_mean)) +
  geom_point() +
  geom_vline(xintercept = 1:42+0.5, color='lightgray') +
  geom_hline(yintercept = med_meanFD, color='red') +
  geom_hline(yintercept = avg_meanFD, color='red', linetype=2) +
  ylab('Mean FD') + xlab('Subjects') +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(angle = 45, size=6, hjust=1)) +
  ggtitle("Mean (Modified) FD")

q_df1$Subject <- factor(q_df1$Subject, levels=rev(sub_by_FDrate))
plt2 <- ggplot(q_df1, aes(x=Subject, y=FD2_rate)) +
  geom_point() +
  geom_vline(xintercept = 1:42+0.5, color='lightgray') +
  geom_hline(yintercept = med_ScrubRate, color='red') +
  geom_hline(yintercept = avg_ScrubRate, color='red', linetype=2) +
  ylab('Proportion of Volumes with FD > 0.2mm') + xlab('Subjects') +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(angle = 45, size=6, hjust=1)) +
  ggtitle("Average Scrubbing Rate (FD > 0.2mm)")

gridExtra::grid.arrange(plt1, plt2, nrow=1)
dev.off()

## Based on the plot above, will use avg_ScrubRate (0.154) as the threshold for high/low movers.
## For the full HCP, we can probably just use 0.15, but should visualize the distribution

# Data remaining after scrubbing -----
# data collected vs data after scrubbing
# for each scrubbing method
# low vs high movers separately visit 1: LR --> RL
# visit 2: RL --> LR (for most subjects)
dcas_df_agg <- NULL
for (tt in seq(length(nT_seq)*2)) {
  use_B <- tt > length(nT_seq) #use second session (RL)?
  tt_B <- if (use_B) { tt - length(nT_seq) } else { NULL }
  nT_tt <- if (use_B) { nT_seq[tt_B] } else { nT_seq[tt] }
  # (42) subjects by (nT_tt) timepoints matrix of FD values from LR1+LR2 concat.
  agFD_tt <- if (use_B) {
    aggFD[idx_RL1,seq(nT_seq[tt_B])] #first X volumes of RL
  } else {
    aggFD[idx_LR1,seq(nT_seq[tt])] #first X volumes of LR
  }

  dcas_df <- FDstats
  # Note which run
  dcas_df$run <- if(use_B) { 'RL' } else { 'LR' }
  # Record the flagging rates
  dcas_df$FD2 <- (nT_tt - rowSums(agFD_tt > .2))/nT_tt
  dcas_df$FD3 <- (nT_tt - rowSums(agFD_tt > .3))/nT_tt
  dcas_df$FD4 <- (nT_tt - rowSums(agFD_tt > .4))/nT_tt
  dcas_df$FD5 <- (nT_tt - rowSums(agFD_tt > .5))/nT_tt
  # Aggregate flagging rates by motion group
  dcas_df <- aggregate(
    dcas_df[paste0("FD", FD_cuts_int)],
    list(dcas_df$motion, dcas_df$run), mean
  )
  colnames(dcas_df)[1:2] <- c("motion","run")
  dcas_df$nT = nT_tt
  dcas_df_agg <- rbind(dcas_df_agg, dcas_df)
}
dcas_df_agg <- tidyr::pivot_longer(dcas_df_agg, seq(3,6))
dcas_df_agg$name <- factor(dcas_df_agg$name, levels=paste0("FD", FD_cuts_int))
dcas_df_agg$run <- factor(dcas_df_agg$run, levels=c('LR','RL'))

colnames(dcas_df_agg)[colnames(dcas_df_agg)=="name"] <- "Scrubbing threshold"
pdf(file.path(dir_plots, "DataCollectedVsRemaining.pdf"), width=7, height=4)
ggplot(
  subset(dcas_df_agg, motion!="Drop"),
  aes(x=nT*0.72/60, y=value, group=`Scrubbing threshold`, color=`Scrubbing threshold`)) +
  geom_hline(yintercept=1) +
  geom_line() + geom_point() +
  scale_color_manual(values=color_gg) +
  xlab("Minutes collected") +
  ylab("Proportion of volumes remaining\n(mean over subjects)") +
  ggtitle("Data collected vs remaining after scrubbing (first visit)") +
  scale_x_continuous(limits=c(2,15), expand=c(0,1), breaks=c(2,5,10,15)) +
  scale_y_continuous(limits=c(.7, 1)) +
  ggthemes::theme_few() + facet_grid(.~motion + run, scales='free_x') + theme(legend.position="bottom")
dev.off()

saveRDS(dcas_df_agg, file = file.path(dir_slate, "results/6_scrubRate/AggFD_byDuration.rds"))

# OLD -----

# [OLD code]
# slen <- dimnames(qt)[[3]]
# pdf(file.path(dir_plots, "ScrubbingRate.pdf"), width=9, height=5)
# for (ss in seq(length(slen))) {
#   q_gg <- as.data.frame.table(arr_gtr[,sub_order,slen[ss],,], responseName="scrub_pct")
#   q_gg$visit <- factor(ifelse(grepl("retest", q_gg$visit), "Retest", "Test"), levels=c("Test", "Retest"))
#   levels(q_gg$FD_cutoff) <- gsub(".*_", "", levels(q_gg$FD_cutoff))
#
#   plt <- ggplot(q_gg, aes(x=subject, y=scrub_pct*100)) +
#     geom_point(aes(color=FD_cutoff, shape=FD_type)) +
#     geom_hline(yintercept=50, lty=2) +
#     geom_hline(yintercept=pct_need, lty=2) +
#     geom_vline(xintercept = 1:42+0.5, color='lightgray') +
#     scale_color_manual(values = c("#e6ab02", color_FD)) +
#     scale_shape_manual(values = c(16, 4)) +
#     ylab('Percent of Volumes Scrubbed') + xlab('Subjects') +
#     scale_y_continuous(breaks=seq(0,100,10)) +
#     ggthemes::theme_few() + facet_grid(.~visit) +
#     scale_y_continuous(expand=c(0,5)) +
#     theme(axis.text.x = element_text(angle = 45, size=6, hjust=1)) +
#     ggtitle(paste0("Scrubbing rates for ", gsub("_" ," ", slen[ss])))
#   print(plt)
# }
# dev.off()

# slen <- dimnames(qt)[[3]]
# pdf(file.path(dir_plots, "ScrubbingRate_mean.pdf"), width=9, height=5)
# for (ss in seq(length(slen))) {
#
#   # TO DO
#   #arr_gt <- lapply(visits_gt, function(z){ colSums(q[z,,"first_14_mins","FD___og_nfc_l4___0.2","FDplus"]) })
#
#   q_gg <- abind::abind(
#     test=colMeans(qty[!grepl("retest", dimnames(qty)[[1]]),sub_order,,,]),
#     retest=colMeans(qty[grepl("retest", dimnames(qty)[[1]]),sub_order,,,]),
#     along=5, use.first.dimnames = TRUE
#   )
#   q_gg <- aperm(q_gg, c(5,1,2,3,4))
#   names(dimnames(q_gg)) <- names(dimnames(qty))
#   q_gg[is.nan(q_gg)] <- 1 # why?
#   q_gg <- as.data.frame.table(q_gg[,,slen[ss],,], responseName="scrub_pct")
#   q_gg$visit <- factor(ifelse(grepl("retest", q_gg$visit), "Retest", "Test"), levels=c("Test", "Retest"))
#   levels(q_gg$FD_cutoff) <- gsub(".*_", "", levels(q_gg$FD_cutoff))
#
#   plt <- ggplot(q_gg, aes(x=subject, y=scrub_pct*100)) +
#     geom_point(aes(color=FD_cutoff, shape=FD_type)) +
#     geom_hline(yintercept=50, lty=2) + geom_vline(xintercept = 1:42+0.5, color='lightgray') +
#     scale_color_manual(values = color_FD) +
#     scale_shape_manual(values = c(16, 4)) +
#     ylab('Percent of Volumes Scrubbed') + xlab('Subjects') +
#     scale_y_continuous(breaks=seq(0,100,10)) +
#     ggthemes::theme_few() + facet_grid(.~visit) +
#     scale_y_continuous(expand=c(0,5)) +
#     theme(axis.text.x = element_text(angle = 45, size=6, hjust=1)) +
#     ggtitle(paste0("Scrubbing rates for ", gsub("_" ," ", slen[ss])))
#   print(plt)
# }
# dev.off()
