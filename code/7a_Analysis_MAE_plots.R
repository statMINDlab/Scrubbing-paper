####################################################################
# Who: Amanda Mejia
# When: Fall 2024-2025
# What: Plot scrubbing rate, MAE, change in MAE, FC bias, and baseline noise
# How: First 7a_Analysis_MAE_compute.R must be run
# Where: Can be run on HPC or locally, since it only requires access to Github files
####################################################################

## ---------------------------------------------------------------------------
source("code/0_SharedCode.R")

library(dplyr)
library(tidyr)
library(ggthemes)
library(cowplot)
library(ggpubr) #geom_bracket
library(patchwork)
library(ggpattern)

#toggle for re-making plots
networklabels <- FALSE

gt <- "FD2"

load(file = file.path(dir_github,'data','subjects_retest.RData')) #sub_drop, sub_himo, sub_lomo (from 6_ScrubbingRate.R)

### HELPER FUNCTIONS

#convert durations from character to numeric
format_duration <- function(x){ as.numeric(gsub("_mins", "", gsub("first_","",x)))*2 }

#convert durations from numeric to "T = X min"
format_duration2 <- function(x){
  x2 <- paste0("T = ", x, " min")
  x2[x2 == "T = 28.44 min"] <- "T = 30 min"
  x2 <- factor(x2, levels = paste0("T = ", seq(5,30,2.5)," min"))
  return(x2)
}

#convert subject IDs to group labels
group_labels <- c('Below-Average Motion', 'Above-Average Motion')
subject_group <- function(IDs, sub_himo){ factor((IDs %in% sub_himo), levels = c(FALSE, TRUE), labels = group_labels) }

#format FD levels
format_FD <- function(x){ factor(x, levels = FD_levels, labels = FD_levels_long) }

####################################################################
### PLOT CENSORING RATE (NOMINAL)
####################################################################

#compute scrubbing rate from flags
flags <- readRDS(file = file.path(dir_github,'results','5_flags','flags.rds'))
flags <- flags[,,,c("None", "FD5", "FD2", "FD2+")] #remove FD4 and FD3
flags_avg <- apply(flags, c(1,2,4), mean, na.rm=TRUE) #average over time points -- proportion of vol's scrubbed
flags_avg <- apply(flags_avg, 2:3, mean, na.rm=TRUE) #average over sessions
names(dimnames(flags_avg)) <- c('Subject', 'FD')
flags_df <- as.data.frame.table(flags_avg, responseName = 'numFlag')
flags_df$FD <- format_FD(flags_df$FD)

#group subjects by motion level
flags_df$group <- subject_group(flags_df$Subject, sub_himo)

#compute means
flags_df_mean <- flags_df %>%
  filter(FD != 'None') %>%
  group_by(group, FD) %>%
  summarize(
    n(),
    mean = mean(numFlag),
    mean_label = paste0(round(100*mean),'%'))

pdf(file.path(dir_github, 'plots', 'scrubRate.pdf'), height=5, width=3.5)
print(ggplot() +
        geom_line(data = subset(flags_df, FD != "None"), aes(FD, numFlag, group = Subject, linetype = group), color = 'gray80', linewidth = 0.5, show.legend = FALSE) +
        #geom_point(data = flags_df_highmotion, aes(x = FD, y = numFlag, group = subject, color = FD), size = 1) +
        geom_line(data = flags_df_mean, aes(x = FD, y = mean, group = group, linetype = group), color = 'black', linewidth = 0.8) +
        geom_point(data = flags_df_mean, aes(x = FD, y = mean, color = FD), size = 1.5) +
        geom_text(data = flags_df_mean, aes(x = FD, y = mean, label = mean_label), vjust = -0.6, hjust = 1.2, size = 2.5) +
        scale_linetype_manual(name = NULL, values = c("Above-Average Motion" = "solid", "Below-Average Motion" = "dashed")) +
        scale_color_manual(values = myCols[-1], guide = 'none') +
        scale_y_continuous(limits = c(0, 0.66), labels = scales::percent) +
        cowplot::theme_cowplot() +
        theme(legend.position = c(0.02, 1.0),
              legend.justification = c("left","top"),
              legend.key.width = unit(1.2, "cm"),
              legend.text = element_text(size = 12),
              legend.key.height = unit(0.5, "cm"),
              axis.text.x=element_text(angle=45,hjust=1),
              axis.title.x = element_blank()) +
        ylab("% Volumes Censored") +
        ggtitle('Censoring Rate'))
dev.off()

p_mae_duration <- setNames(vector(mode='list', length=nB), baseNames)
p_dmae_duration <- setNames(vector(mode='list', length=nB), baseNames)
p_dmae_durbymot <- setNames(vector(mode='list', length=nB), baseNames)
p_blSD_bpeff <- setNames(vector(mode='list', length=nB), baseNames)

# Loop over P36, FIX, FIX_GSR
for (bb in 2:nB) {

  base <- baseNames[bb]

  print(paste0('~~~~~~~~~~~~~~~~~~~~~~~~~~ ',base, ' ~~~~~~~~~~~~~~~~~~~~~~~~~~ '))

  ####################################################################
  #### PLOT EFFECTIVE CENSORING RATE
  ####################################################################

  load(file = file.path(dir_github, 'results','7_MAE',paste0('T_eff_sum_',base,'.RData'))) #T_eff_sum

  #group by motion
  T_eff_sum <- filter(T_eff_sum, FD %in% FD_levels)
  T_eff_sum$FD <- format_FD(T_eff_sum$FD)
  T_eff_sum$group <- subject_group(T_eff_sum$subject, sub_himo)

  #compute change in effective scan duration
  T_eff_sum0 <- T_eff_sum %>% filter(FD == 'None') %>% select(subject, nVol_eff0 = nVol_eff)
  T_eff_sum <- T_eff_sum %>% left_join(T_eff_sum0, by = 'subject') %>%
    mutate(loss_Teff = nVol_eff0 - nVol_eff,
           pct_loss_Teff = (nVol_eff0 - nVol_eff) / nVol_eff0)

  #compute mean over subjects
  T_eff_sum_mean <- T_eff_sum %>%
    filter(FD != 'None') %>%
    group_by(group, FD) %>%
    summarize(mean_loss = mean(loss_Teff),
              mean_pct_loss = mean(pct_loss_Teff),
              mean_pct_loss_label = paste0(round(mean(pct_loss_Teff)*100),'%'))

  #make plot
  pdf(file.path(dir_github, 'plots', paste0("scrubRateEff_", base, '.pdf')), height=4, width=3.6)
  print(ggplot() +
          geom_line(data = subset(T_eff_sum, FD != "None"), aes(FD, pct_loss_Teff, group = subject, linetype = group), color = 'gray80', linewidth = 0.5, show.legend = FALSE) +
          #geom_point(data = flags_df_highmotion, aes(x = FD, y = numFlag, group = subject, color = FD), size = 1) +
          geom_line(data = T_eff_sum_mean, aes(x = FD, y = mean_pct_loss, group = group, linetype = group), color = 'black', linewidth = 0.8) +
          geom_point(data = T_eff_sum_mean, aes(x = FD, y = mean_pct_loss, color = FD), size = 1.5) +
          geom_text(data = T_eff_sum_mean, aes(x = FD, y = mean_pct_loss, label = mean_pct_loss_label), vjust = -0.4, hjust = 1.2, size = 3) +
          scale_linetype_manual(name = NULL, values = c("Above-Average Motion" = "solid", "Below-Average Motion" = "dashed")) +
          scale_color_manual(values = myCols[-1], guide = 'none') +
          scale_x_discrete(expand=expansion(add=c(.25, .15))) +
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limits = c(0,0.66), labels = scales::percent) +
          cowplot::theme_cowplot() +
          theme(legend.position = c(0.02, 0.95),
                legend.justification = c("left","top"),
                legend.key.width = unit(1.2, "cm"),
                legend.text = element_text(size = 12),
                legend.key.height = unit(0.5, "cm"),
                axis.text.x=element_text(angle=45,hjust=1),
                axis.title.x = element_blank()) +
          #ggtitle('Effective Censoring Rate')) + # Hidden
          ylab("% Reduction")) # % Loss in Effective Scan Duration
  dev.off()

  ####################################################################
  ### PLOT FC ERROR VS SCAN DURATION
  ####################################################################

  ### TWO VERSIONS: including all subjects and split into high/low motion
  ### DURATION: Vary duration for T=5 to T=30 min

  load(file = file.path(dir_github, 'results','7_MAE',paste0('MAE_',base,'.RData'))) #MAE_subj, MAE_edge_summ
  MAE_edge_summ$FD <- format_FD(MAE_edge_summ$FD)
  MAE_edge_summ$Duration <- format_duration(MAE_edge_summ$Duration)
  MAE_edge_summ$group <- factor(MAE_edge_summ$group, levels = c('all','lomo', 'himo'), labels = c('all',group_labels))

  #ungrouped
  pdf(file.path(dir_github, "plots", "MAE", base, "MAE_duration.pdf"), height=4, width=4.5)
  p_mae_duration[[bb]] <- ggplot(filter(MAE_edge_summ, group=='all'), aes(x=Duration, y=sqrt(MSE))) +
      geom_line(aes(group=FD, color=FD)) + geom_point(aes(color=FD)) + # removed from geom_point: linewidth=2 (not a valid param for geom_point)
      scale_color_manual('Censoring\nMethod', values = myCols) +
      scale_x_continuous(breaks=c(5,10,20,30), limits=c(3,30)) +
      xlab("Scan Duration (min)") +
      ylab("Root Mean Squared Error (rMSE)") +
      ggtitle("FC Error") +
      cowplot::theme_cowplot() +
      theme(legend.title = element_blank(),
            legend.position = c(0.95, 0.95),
            legend.justification = c("right", "top"))
  print(p_mae_duration[[bb]])
  dev.off()

  #grouped
  pdf(file.path(dir_github, "plots", "MAE", base, "MAE_duration_bymotion.pdf"), height=4, width=5)
  print(ggplot(filter(MAE_edge_summ, group != 'all'), aes(x=Duration, y=sqrt(MSE))) +
          geom_line(aes(group=FD, color=FD)) + geom_point(aes(color=FD)) + # removed from geom_point: linewidth=2
          scale_color_manual('Censoring\nMethod', values = myCols) +
          scale_x_continuous(breaks=c(5,10,20,30), limits=c(3,30)) +
          xlab("Scan Duration (min)") +
          ylab("Root Mean Squared Error") +
          ggtitle("FC Error") +
          cowplot::theme_cowplot() +
          facet_grid(. ~ group) +
          theme(legend.title = element_blank(),
                legend.position = c(0.3, 0.95),
                legend.justification = c("right", "top"),
                strip.background = element_blank(),
                strip.text = element_text(face = 'bold')))
  dev.off()

  #by subject ("boxplots")
  MAE_subj$Duration <- format_duration(MAE_subj$Duration)
  MAE_subj$Duration2 <- format_duration2(MAE_subj$Duration)
  MAE_subj$FD <- format_FD(MAE_subj$FD)
  MAE_subj$group <- subject_group(MAE_subj$Subject, sub_himo)
  df_plot <- subset(MAE_subj, Duration %in% c(5, 10, 15, 20, 25, 28.44))
  pdf(file.path(dir_github, "plots", "MAE", base, "MAE_boxplots.pdf"), height=5, width=10)
  print(ggplot(df_plot) +
          geom_line(data = df_plot, aes(x=FD, y=rMSE, group=Subject), alpha=0.2) +
          geom_point(aes(x=FD, y=rMSE, group=FD, color=FD, shape=group), alpha=0.7, size=2) +
          facet_grid(. ~ Duration2) +
          scale_color_manual(values = myCols, guide="none") +
          xlab("Censoring Level") +
          ylab("Root Mean Squared Error") +
          ggtitle("FC Error") +
          cowplot::theme_cowplot() +
          theme(legend.title=element_blank(),
                legend.position='bottom',
                legend.justification="center",
                axis.text.x=element_text(angle=45,hjust=1),
                strip.background = element_blank(),
                strip.text = element_text(face = 'bold')))
  dev.off()


  ####################################################################
  ### PLOT CHANGE IN FC ERROR VS DURATION
  ####################################################################

  MAE_diff_subj <- readRDS(file = file.path(dir_github, 'results','7_MAE',paste0('MAE_diff_subj_',base,'.rds')))

  ### PLOT 1: By subject ("boxplots") with significance testing

  #MAE difference by subject
  MAE_diff_subj$Duration <- format_duration(MAE_diff_subj$Duration)
  MAE_diff_subj$Duration2 <- format_duration2(MAE_diff_subj$Duration)
  MAE_diff_subj$FD <- format_FD(MAE_diff_subj$FD)
  MAE_diff_subj$group <- subject_group(MAE_diff_subj$Subject, sub_himo)
  df_plot <- subset(MAE_diff_subj, (FD != 'None') & (Duration %in% c(5, 10, 15, 20, 25, 28.44)))

  #perform Wilcox tests
  durations <- c(5, 10, 15, 20, 25, 28.44)
  pvals <- matrix(NA, nrow=length(durations), ncol=3)
  rownames(pvals) <- durations
  col_names <- c()
  which_col <- 1
  #loop over censoring levels
  for(s1 in 1:3){
    s2 <- s1+1
    #compare censoring levels
    FD_s1 <- FD_levels_long[s1]
    FD_s2 <- FD_levels_long[s2]
    col_names <- c(col_names, paste0(FD_s1,"-",FD_s2))
    for(dd in 1:length(durations)){
      dur <- durations[dd]
      MAE_diff_subj_sd1 <- filter(MAE_diff_subj, Duration==dur, FD==FD_s1)
      MAE_diff_subj_sd2 <- filter(MAE_diff_subj, Duration==dur, FD==FD_s2)
      pvals[dd,which_col] <- wilcox.test(MAE_diff_subj_sd1$rMSE_delta,
                                         MAE_diff_subj_sd2$rMSE_delta, paired = TRUE)$p.value
    }
    which_col <- which_col + 1
  }
  colnames(pvals) <- col_names
  pvals_adj <- pvals*length(pvals) #bonferroni correction
  names(dimnames(pvals_adj)) <- c("Duration","Comparison")
  pvals_sig_df <- as.data.frame.table(pvals_adj, responseName = "pvalue")
  pvals_sig_df$Duration2 <- format_duration2(pvals_sig_df$Duration)
  pvals_sig_df <- filter(pvals_sig_df, pvalue < 0.05)
  pvals_sig_df$Star <- ifelse(pvals_sig_df$pvalue < 0.001, '***', ifelse(pvals_sig_df$pvalue < 0.01, '**', '*'))
  pvals_sig_df0 <- filter(pvals_sig_df, Comparison=="None-Lenient") #not shown on plot
  pvals_sig_df1 <- filter(pvals_sig_df, Comparison=="Lenient-Stringent")
  pvals_sig_df2 <- filter(pvals_sig_df, Comparison=="Stringent-Expanded")

  print(pvals_sig_df0) #to report results of Lenient vs None in text

  #height of significance bars
  ymin <- min(0,min(df_plot$rMSE_delta)) #limit of plot
  ymax <- max(df_plot$rMSE_delta)*1.2 #limit of plot
  y1 <- ymax*0.85
  y2 <- ymax*0.95

  #PLOT 1, VERSION A
  pdf(file.path(dir_github, "plots", "MAE", base, "dMAE_boxplots_new.pdf"), height=5, width=10)
  p <- (ggplot(df_plot) +
          geom_line(aes(x=FD, y=rMSE_delta, group=Subject), alpha=0.2) +
          geom_point(aes(x=FD, y=rMSE_delta, group=FD, color=FD, shape=group), alpha=0.7, size=2) +
          geom_hline(yintercept = 0, linetype=2) +
          scale_shape_manual('Motion Level', values=c(19,8)) +
          facet_grid(. ~ Duration2) +
          scale_color_manual(values = myCols[-1], guide="none") +
          scale_y_continuous(labels = scales::percent, limits=c(ymin,ymax)) +
          xlab("Censoring Level") + ylab("Change in rMSE") +
          ggtitle("Change in FC Error") +
          cowplot::theme_cowplot() +
          theme(legend.title=element_blank(),
                legend.position='bottom',
                legend.justification="center",
                axis.text.x=element_text(angle=45,hjust=1),
                strip.background = element_blank(),
                strip.text = element_text(face = 'bold')))
  print(p + geom_bracket(data = pvals_sig_df1, aes(xmin=1, xmax=2, y.position=y1, label=Star), inherit.aes = FALSE) + #Lenient-Stringent
          geom_bracket(data = pvals_sig_df2, aes(xmin=2, xmax=3, y.position=y2, label=Star), inherit.aes = FALSE)) #Stringent-Expanded
  dev.off()

  #PLOT 1, VERSION B (just T = 5 min and T = 30 min)
  df_plot <- subset(df_plot, (Duration %in% c(5, 28.44)))
  pvals_sig_df1 <- subset(pvals_sig_df1, (Duration %in% c(5, 28.44)))
  pvals_sig_df2 <- subset(pvals_sig_df2, (Duration %in% c(5, 28.44)))
  pdf(file.path(dir_github, "plots", "MAE", base, "dMAE_boxplots_T5T30.pdf"), height=5, width=7)
  p <- (ggplot(df_plot) +
          geom_line(aes(x=FD, y=rMSE_delta, group=Subject), alpha=0.2) +
          geom_point(aes(x=FD, y=rMSE_delta, group=FD, color=FD, shape=group), alpha=0.7, size=2) +
          geom_hline(yintercept = 0, linetype=2) +
          scale_shape_manual('Motion Level', values=c(19,8)) +
          facet_grid(. ~ Duration2) +
          scale_color_manual(values = myCols[-1], guide="none") +
          scale_y_continuous(labels = scales::percent, limits=c(ymin,ymax)) +
          xlab("Censoring Level") + ylab("Change in rMSE") +
          ggtitle("Change in FC Error") +
          cowplot::theme_cowplot() +
          theme(legend.title=element_blank(),
                legend.position='bottom',
                legend.justification="center",
                axis.text.x=element_text(angle=45,hjust=1),
                strip.background = element_blank(),
                strip.text = element_text(face = 'bold')))
  print(p + geom_bracket(data = pvals_sig_df1, aes(xmin=1, xmax=2, y.position=y1, label=Star), inherit.aes = FALSE) + #Lenient-Stringent
          geom_bracket(data = pvals_sig_df2, aes(xmin=2, xmax=3, y.position=y2, label=Star), inherit.aes = FALSE)) #Stringent-Expanded
  dev.off()


  ### PLOT 2: Mean change in rMSE over subjects

  #ungrouped
  MAE_diff_summ <- MAE_diff_subj %>% filter(FD != 'None') %>% group_by(Duration, FD) %>% summarize(rMSE_delta = mean(rMSE_delta))
  ymax <- max(MAE_diff_summ$rMSE_delta)*1.25
  pdf(file.path(dir_github, "plots", "MAE", base, 'dMAE_duration.pdf'), height=4, width=4.5)
  p_dmae_duration[[bb]] <- ggplot(MAE_diff_summ, aes(x=Duration, y=rMSE_delta)) +
          geom_line(aes(group=FD, color=FD)) + geom_point(aes(color=FD)) + # removed from geom_point: linewidth=2
          scale_color_manual('Censoring\nMethod', values = myCols[-1]) +
          scale_x_continuous(breaks=c(5,10,20,30), limits=c(3,30)) +
          scale_y_continuous(labels = scales::percent, limits=c(0,ymax)) +
          geom_hline(yintercept = 0, linetype=2) +
          xlab("Scan Duration (min)") +
          ylab("Change in rMSE") +
          ggtitle("Change in FC Error") +
          cowplot::theme_cowplot() +
          theme(legend.title = element_blank(),
                legend.position = c(0.95, 0.99),
                legend.justification = c("right", "top"))
  print(p_dmae_duration[[bb]])
  dev.off()

  #grouped by motion level
  MAE_diff_summ_grouped <- MAE_diff_subj %>% filter(FD != 'None') %>% group_by(Duration, FD, group) %>% summarize(rMSE_delta = mean(rMSE_delta))
  ymax <- max(MAE_diff_summ_grouped$rMSE_delta)*1.25
  pdf(file.path(dir_github, "plots", "MAE", base, 'dMAE_duration_bymotion.pdf'), height=4, width=5)
  p_dmae_durbymot[[bb]] <- ggplot(MAE_diff_summ_grouped, aes(x=Duration, y=rMSE_delta)) +
          geom_line(aes(group=FD, color=FD)) + geom_point(aes(color=FD)) + # rm: linewidth=2
          scale_color_manual('Censoring\nMethod', values = myCols[-1]) +
          scale_x_continuous(breaks=c(5,10,20,30), limits=c(3,30)) +
          scale_y_continuous(labels = scales::percent) +
          geom_hline(yintercept = 0, linetype=2) +
          xlab("Scan Duration (min)") +
          ylab("Change in rMSE") +
          ggtitle("Change in FC Error") +
          cowplot::theme_cowplot() +
          facet_grid(. ~ group, labeller = labeller(group = as_labeller(group_labels))) +
          theme(legend.title = element_blank(),
                legend.position = c(0.3, 0.95),
                legend.justification = c("right", "top"),
                strip.background = element_blank(),
                strip.text = element_text(face = 'bold'))
  plot(p_dmae_durbymot[[bb]])
  dev.off()


  ####################################################################
  ### PLOT DURATION VS. MAE (reverse of above, "budget inflation" analysis)
  ####################################################################

  ### Plot a few subjects/edges to illustrate approach

  #shows that FD5 itself can sometimes achieve the target MAE with shorter duration
  #so we use that duration to compute the % change in duration for other scrubbing levels
  #also shows that the required scan duration can be > 30 minutes

  ## by Edges

  load(file = file.path(dir_github, 'results','7_MAE',paste0('minDurationIllustration_',base,'.RData'))) #edges_illustration, MAE_edge_illustration, MAE_minDur_edge_illustration
  MAE_edge_illustration$FD <- format_FD(MAE_edge_illustration$FD)
  MAE_minDur_edge_illustration$FD <- format_FD(MAE_minDur_edge_illustration$FD)

  label_both2 <- function(labels){
    label_both(labels, multi_line=TRUE, sep=' ')
  }

  pdf(file.path(dir_github, "plots", "MAE", base, 'durationChange_illustration_edges.pdf'), height = 7, width = 5)
  annotate_loc <- MAE_edge_illustration %>%
    filter(Edge == edges_illustration[1]) %>% group_by(Edge) %>%
    summarize(ymax = max(rMSE), FD = 'Lenient') #where to locate "legend" for blue dot
  print(ggplot(MAE_edge_illustration, aes(x=Duration, y=rMSE, color=FD, group=FD)) +
          geom_point() + geom_line(alpha = 0.5) +
          geom_vline(data=MAE_minDur_edge_illustration, aes(xintercept = minDur, group=FD, color=FD), lty=2, alpha=0.5) + #minDur for each FD
          geom_hline(data = MAE_minDur_edge_illustration, aes(yintercept = target)) +
          geom_point(data = filter(MAE_minDur_edge_illustration, FD=='Lenient'), x = 17.5, aes(y = target), color='black', fill='deepskyblue', size = 3, shape = 21) +
          geom_point(data = annotate_loc, aes(y = ymax-0.001), x=24, color='black', fill='deepskyblue', size = 3, shape = 21) +
          geom_text(data = annotate_loc, aes(y = ymax-0.001), x=27, label = 'Target rMSE', color='black', size = 3) +
          scale_color_manual(name='', values = myCols[-1]) + ylab('Root Mean Squared Error') +
          facet_grid(Edge ~ ., scales = 'free', labeller = label_both2) +
          theme_few() + theme(legend.position = 'bottom'))
  dev.off()

  ## by Subjects

  load(file = file.path(dir_github, 'results','7_MAE',paste0('minDuration_subj_',base,'.RData'))) #MAE_minDur_subj
  MAE_minDur_subj$FD <- format_FD(MAE_minDur_subj$FD)
  MAE_minDur_subj$group <- subject_group(MAE_minDur_subj$Subject, sub_himo)

  set.seed(974989)
  subjects_illustration <- sort(sample(unique(MAE_subj$Subject), 3))
  MAE_minDur_subj_illustration <- filter(MAE_minDur_subj, Subject %in% subjects_illustration) #for vertical lines indicating required durations
  MAE_subj_illustration <- filter(MAE_subj, Subject %in% subjects_illustration, FD != 'None') #for rMSE lines
  annotate_loc <- MAE_subj_illustration %>%
    filter(Subject == subjects_illustration[1]) %>% group_by(Subject) %>%
    summarize(ymax = max(rMSE), FD = 'Lenient') #where to locate "legend" for blue dot

  pdf(file.path(dir_github, "plots", "MAE", base, 'durationChange_illustration_subjects.pdf'), height = 7, width = 5)
  print(ggplot(MAE_subj_illustration, aes(x=Duration, y=rMSE, color=FD, group=FD)) +
          geom_point() + geom_line(alpha = 0.5) +
          geom_vline(data= MAE_minDur_subj_illustration, aes(xintercept = minDur, group=FD, color=FD), lty=2, alpha=0.5) + #minDur for each FD
          geom_hline(data = MAE_minDur_subj_illustration, aes(yintercept = target)) +
          geom_point(data = filter(MAE_minDur_subj_illustration, FD=='Lenient'), x = 17.5, aes(y = target), color='black', fill='deepskyblue', size = 3, shape = 21) +
          geom_point(data = annotate_loc, aes(y = ymax-0.001), x=24, color='black', fill='deepskyblue', size = 3, shape = 21) +
          geom_text(data = annotate_loc, aes(y = ymax-0.001), x=27, label = 'Target rMSE', color='black', size = 3) +
          scale_color_manual(name = '', values = myCols[-1]) + ylab('Root Mean Squared Error') +
          facet_grid(Subject ~ ., scales = 'free', labeller = label_both2) +
          theme_few() + theme(legend.position = 'bottom'))
  dev.off()

  ### Plot change in required scan duration by subject

  #compute means
  MAE_minDur_mean <- MAE_minDur_subj %>%
    filter(FD != 'Lenient') %>%
    group_by(group, FD) %>%
    summarize(
      mean_mins = mean(Duration_change_mins),
      mean_mins_label = round(mean_mins, 1),
      mean_pct = mean(Duration_change_pct),
      mean_pct_label = paste0(round(100*mean_pct), '%'))

  pdf(file.path(dir_github, "plots", "MAE", base, 'durationChange_subjects.pdf'), height = 6, width = 4.5)
  p <- ggplot(filter(MAE_minDur_subj, FD != 'Lenient'), aes(x = FD, y = Duration_change_mins)) +
    geom_line(aes(group = Subject, linetype = group), alpha=0.1) +
    geom_line(data = MAE_minDur_mean, aes(x = FD, y = mean_mins, group = group, linetype = group), color = 'black', linewidth = 0.8) +
    geom_point(data = MAE_minDur_mean, aes(x = FD, y = mean_mins, color = FD), size = 3) +
    geom_text(data = MAE_minDur_mean, aes(x = FD, y = mean_mins, label = mean_mins_label), vjust = -0.7, hjust = 1.3, size = 4) +
    scale_color_manual(values = myCols[-(1:2)], guide='none') +
    scale_linetype_manual(values = c(2,1)) +
    scale_x_discrete(expand = c(0.15,0.15)) +
    xlab("Censoring Level") +
    ylab("Change in Required Scan Duration (mins)") +
    ggtitle("Scan Duration to Maintain FC \nAccuracy (vs Lenient Censoring)") +
    cowplot::theme_cowplot() +
    theme(legend.title=element_blank(),
          legend.position = c(0.05, 0.95),
          legend.justification = c("left","top"),
          legend.key.width = unit(2, "line"))
  if(base != 'FIX') p <- p + ylim(-0.1,20)
  print(p)
  p <- ggplot(filter(MAE_minDur_subj, FD != 'Lenient'), aes(x = FD, y = Duration_change_pct)) +
          geom_line(aes(group = Subject, linetype = group), alpha=0.1) +
          geom_line(data = MAE_minDur_mean, aes(x = FD, y = mean_pct, group = group, linetype = group), color = 'black', linewidth = 0.8) +
          geom_point(data = MAE_minDur_mean, aes(x = FD, y = mean_pct, color = FD), size = 3) +
          geom_text(data = MAE_minDur_mean, aes(x = FD, y = mean_pct, label = mean_pct_label), vjust = -0.7, hjust = 1.1, size = 4) +
          scale_color_manual(values = myCols[-(1:2)], guide='none') +
          scale_linetype_manual(values = c(2,1)) +
          scale_y_continuous(labels = scales::percent) +
          scale_x_discrete(expand = c(0.15,0.15)) +
          xlab("Censoring Level") +
          ylab("% Change in Required Scan Duration") +
          ggtitle("Scan Duration to Maintain FC \nAccuracy (vs Lenient Censoring)") +
          cowplot::theme_cowplot() +
          theme(legend.title=element_blank(),
                legend.position = c(0.05, 0.95),
                legend.justification = c("left","top"),
                legend.key.width = unit(2, "line"))
  if(base != 'FIX') p <- p + scale_y_continuous(labels = scales::percent, limits = c(-0.01,1.7))
  print(p)
  dev.off()


  ####################################################################
  ### PLOTS OF FC BIAS
  ####################################################################

  load(file = file.path(dir_github, 'results','7_MAE',paste0('bias_',base,'.RData'))) #bias_subj

  bias_subj$FD <- format_FD(bias_subj$FD)
  bias_subj$Duration <- format_duration(bias_subj$Duration)
  bias_subj$group <- subject_group(bias_subj$Subject, sub_himo)

  bias_summ <- bias_subj %>% group_by(Duration, FD, group) %>% summarize(bias_mean = mean(bias))
  pdf(file.path(dir_github, "plots", "MAE", base, 'bias_duration.pdf'), height=4, width=5)
  print(ggplot(bias_summ, aes(x=Duration, y=bias_mean)) +
          geom_line(aes(group=FD, color=FD)) + geom_point(aes(color=FD)) + # linewidth=2
          scale_color_manual('Censoring\nMethod', values = myCols) +
          scale_x_continuous(breaks=c(5,10,20,30), limits=c(3,30)) +
          scale_y_continuous(expand = c(0.02,0.02)) +
          xlab("Scan Duration (min)") +
          ylab("Bias in FC") +
          ggtitle("FC Bias vs. Ground Truth") +
          cowplot::theme_cowplot() +
          facet_grid(. ~ group, labeller = labeller(group = as_labeller(group_labels))) +
          theme(legend.title = element_blank(),
                legend.position = c(0.95, 0.95),
                legend.justification = c("right", "top"),
                strip.background = element_blank(),
                strip.text = element_text(face = 'bold')))
  dev.off()

  ####################################################################
  ### PLOTS OF BASELINE NOISE
  ####################################################################

  load(file = file.path(dir_github, 'results','7_MAE',paste0('baselineSD_',base,'.RData'))) #baselineSD_subj

  plot_baselineSD <- function(df, colors, ylim, title, relabel=FALSE){

    if (relabel) {
      df$group <- factor(
        df$group,
        levels=levels(df$group),
        labels=gsub(" ", "\n", levels(df$group))
      )
    }

    ggplot(filter(df, FD != "None"),
           aes(x = FD, y = baselineSDeff_mean_delta)) +
      #geom_line(aes(group=Subject, linetype=group), alpha=0.1) +
      #geom_boxplot_pattern(aes(group = interaction(group, FD), pattern = group, fill = FD),  pattern_spacing = 0.01, outliers = FALSE) +
      geom_boxplot(aes(group = FD, fill = FD), outliers = FALSE) +
      scale_fill_manual(values = colors, guide = 'none') +
      #scale_pattern_manual(values = c('none','stripe')) +
      #geom_point(size=2) +
      ylab("% Change") + # Change in Noise vs. No Censoring
      ggtitle(title) +
      cowplot::theme_cowplot() +
      scale_y_continuous(labels = scales::percent) +
      coord_cartesian(ylim = ylim) +
      facet_grid(. ~ group) +
      theme(legend.position = 'none',
            plot.title = element_blank(), # for Figure
            # plot.title = element_text(hjust = 0.5), # for Figure
            strip.text = element_text(face = "plain"), # for Figure
            # strip.text = element_text(face = 'bold') # for Figure
            #legend.position = c(0.05,0.1),
            #legend.title = element_blank(),
            axis.text.x=element_text(angle=45,hjust=1),
            axis.title.x = element_blank(),
            strip.background = element_blank())
  }

  ### Version 1: All scrubbing levels

  pdf(file.path(dir_github, "plots", "MAE", base, "baselineSD_boxplot_eff.pdf"), height=5, width=5)
  print(plot_baselineSD(df = baselineSD_subj, colors=myCols2[-1], ylim=c(-0.35, 0.15), title = baseNames_long[bb]))
  dev.off()

  ### Version 2: Lenient/Stringent/Expanded

  baselineSD_subj <- filter(baselineSD_subj, FD %in% FD_levels)
  baselineSD_subj$FD <- format_FD(baselineSD_subj$FD)

  pdf(file.path(dir_github, "plots", "MAE", base, "baselineSD_boxplot_eff2.pdf"), height=4.5, width=3.6)
  p_blSD_bpeff[[bb]] <- plot_baselineSD(
    df = baselineSD_subj, colors=myCols[-1], ylim=c(-0.35, 0.15), title = baseNames_long[bb], relabel=TRUE
  )
  print(p_blSD_bpeff)
  dev.off()

  ### Version 3: Line plot for Lenient/Stringent/Expanded

  #compute mean over subjects
  baselineSD_subj_lomo <- baselineSD_subj %>%
    filter(FD != "None", group == "Below-Average Motion", !(base == "FIX_GSR" & Subject %in% c("662551"))) %>%
    group_by(group, FD) %>%
    summarize(mean_delta = mean(baselineSDeff_mean_delta)) %>%
    mutate(mean_delta_label = ifelse(round(mean_delta*100) > 0, paste0("+",round(mean_delta*100),"%"), paste0(round(mean_delta*100),"%")))

  baselineSD_subj_himo <- baselineSD_subj %>%
    filter(FD != "None", group == "Above-Average Motion", !(base == "FIX_GSR" & Subject %in% c("662551"))) %>%
    group_by(group, FD) %>%
    summarize(mean_delta = mean(baselineSDeff_mean_delta)) %>%
    mutate(mean_delta_label = ifelse(round(mean_delta*100) > 0, paste0("+",round(mean_delta*100),"%"), paste0(round(mean_delta*100),"%")))

  top_labels <- if (base == "P36") { baselineSD_subj_lomo } else { baselineSD_subj_himo }
  bottom_labels <- if (base == "P36") { baselineSD_subj_himo } else { baselineSD_subj_lomo }
  top_vnudge <- if (base == "FIX_GSR") { -0.1 } else { -1.0 }
  bottom_vnudge <- if (base == "P36") { 1.0 } else if (base == "FIX") { 1.4 } else { 1.6 }
  top_hnudge <- if (base == "FIX_GSR") { -0.3 } else { 0.5 }
  bottom_hnudge <- if (base == "FIX_GSR") { 0.5 } else { 1.2 }

  pdf(file.path(dir_github, 'plots', "MAE", base, "baselineSD_lineplot.pdf"), height=4, width=3.6)
  print(ggplot() +
          geom_line(data = subset(baselineSD_subj, FD != "None"), aes(FD, baselineSDeff_mean_delta, group = Subject, linetype = group), color = 'gray80', linewidth = 0.5, show.legend = FALSE) +
          geom_line(data = baselineSD_subj_lomo, aes(x = FD, y = mean_delta, group = group, linetype = group), color = 'black', linewidth = 0.8) +
          geom_line(data = baselineSD_subj_himo, aes(x = FD, y = mean_delta, group = group, linetype = group), color = 'black', linewidth = 0.8) +

          geom_point(data = baselineSD_subj_lomo, aes(x = FD, y = mean_delta, color = FD), size = 1.5) +
          geom_point(data = baselineSD_subj_himo, aes(x = FD, y = mean_delta, color = FD), size = 1.5) +

          geom_text(data = top_labels, aes(x = FD, y = mean_delta, label = mean_delta_label), vjust = top_vnudge, hjust = top_hnudge, size = 2.80) +
          geom_text(data = bottom_labels, aes(x = FD, y = mean_delta, label = mean_delta_label), vjust = bottom_vnudge, hjust = bottom_hnudge, size = 2.80) +
          scale_linetype_manual(name = NULL, values = c("Above-Average Motion" = "solid", "Below-Average Motion" = "dashed")) +
          scale_color_manual(values = myCols[-1], guide = 'none') +
          #scale_x_discrete(expand=expansion(add=c(.25, .15))) +
          #scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limits = c(0,0.66), labels = scales::percent) +
          scale_y_continuous(limits = c(-0.35, 0.15), labels = scales::percent) +
          cowplot::theme_cowplot() +
          theme(legend.position = c(0.02, 0.02),
                legend.justification = c("left","bottom"),
                legend.key.width = unit(1.2, "cm"),
                legend.text = element_text(size = 10),
                legend.key.height = unit(0.5, "cm"),
                axis.text.x=element_text(angle=45,hjust=1),
                axis.title.x = element_blank()) +
          ylab("% Change"))
  dev.off()

} #end loop over baseline denoising

# Pre-comps
## Figure 1
pdf(file.path(dir_github, "figures", "precomps", "fig1_AB.pdf"), height=8, width=8.5)
p <- cowplot::plot_grid(
  p_mae_duration$P36 + theme(plot.title = element_blank(), axis.title.y.left = element_text(angle = 90, vjust = 0.5)),
  NULL,
  p_dmae_duration$P36 + theme(plot.title = element_blank(), axis.title.y.left = element_text(angle = 90, vjust = 0.5)),
  p_mae_duration$FIX_GSR + theme(plot.title = element_blank(), axis.title.y.left = element_text(angle = 90, vjust = 0.5)),
  NULL,
  p_dmae_duration$FIX_GSR + theme(plot.title = element_blank(), axis.title.y.left = element_text(angle = 90, vjust = 0.5)),
  align = "v", axis = "n", # align y-axis titles
  rel_widths=c(.8, .1, .8),
  nrow=2
)
print(p)
dev.off()

pdf(file.path(dir_github, "figures", "precomps", "fig1_AB_legend.pdf"), height=3, width=5)
p <- p_mae_duration$P36 +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.direction = "horizontal")
p <- get_legend(p)
grid::grid.draw(p)
dev.off()

# Pre-comps
## Figure 2
pdf(file.path(dir_github, "figures", "precomps", "fig2_A.pdf"), height=4, width=10)
p <- cowplot::plot_grid(
  p_dmae_durbymot$P36 + theme(plot.title = element_blank(), strip.text = element_text(face = "plain")),
  p_dmae_durbymot$FIX_GSR + theme(plot.title = element_blank(), strip.text = element_text(face = "plain")),
  align="hv", axis="lb", nrow=1
)
print(p)
dev.off()

pdf(file.path(dir_github, "figures", "precomps", "fig2_D.pdf"), height=4, width=7)
p <- cowplot::plot_grid(
  p_blSD_bpeff$P36 + theme(
    axis.line = element_line(colour = "black", linewidth = .55),
    axis.ticks = element_line(linewidth = .55)),
  p_blSD_bpeff$FIX_GSR + theme(
    axis.line = element_line(colour = "black", linewidth = .55),
    axis.ticks = element_line(linewidth = .55)),
  align="hv", axis="lb", nrow=1
)
print(p)
dev.off()

