####################################################################
# Who: Amanda Mejia
# When: Fall 2024-2025
# What: Visualize edge-wise MAE, change in MAE, FC bias, and baseline noise
# How: First 7a_Analysis_MAE_compute.R must be run
# Where: Can be run on HPC or locally, since it only requires access to Github files
####################################################################

## ---------------------------------------------------------------------------
source("code/0_SharedCode.R")

library(dplyr)

source(file.path(dir_github, "code", "FC_vis_funs.R")) #tools and functions for FC matrix visualization

networklabels <- FALSE

### HELPER FUNCTIONS

#format FD levels
format_FD <- function(x){ factor(x, FD_levels, labels = FD_levels_long) }

#this function loops over scrubbing levels (typically to make a multi-page PDF)
plot_FC_mat <- function(df, val_name, legTitle, title, FD, zlim = c(-15, 15), colFUN = seq3div){
  #loop over scrubbing (reordered from least to most)
  df <- as.data.frame(df)
  vals_ss <- df[df$FD==FD,val_name]
  vals_ss[vals_ss < zlim[1]] <- zlim[1]
  vals_ss[vals_ss > zlim[2]] <- zlim[2]
  plt <- ggcorrplot2(
    values = vals_ss,
    colFUN = colFUN,
    divColor="black",
    title = paste0(title, ", ", FD),
    legTitle = legTitle,
    lim=zlim)
  plt <- plt + labs(x = "Node", y = "Node") +
      theme_classic() + theme(
        legend.position = "bottom",
        axis.title.x = element_text(size = 24, margin = margin(t = 10)),
        axis.title.y = element_text(size = 24, margin = margin(t = 10)),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
  return(plt)
}

baseNames <- baseNames[-1] # P36, FIX, FIX_GSR
for (base in baseNames) {

  print(paste0('~~~~~~~~~~~~~~~~~~~~~~~~~~ ',base, ' ~~~~~~~~~~~~~~~~~~~~~~~~~~ '))


  ##########################################################################
  ### CHANGE IN RMSE

  MAE_edge <- readRDS(file = file.path(dir_results, 'results','7_MAE',paste0('MAE_diff_edge_',base,'.rds')))
  MAE_edge$rMSE_delta <- MAE_edge$rMSE_delta*100
  MAE_edge$FD <- format_FD(MAE_edge$FD)

  for(ss in FD_levels_long[2:4]){
    fname_ss <- paste0('dMAE_', ss, '.png')
    png(file.path(dir_github, "plots", "MAE", base, fname_ss), width=1000, height=900, res=200)
    print(plot_FC_mat(df = filter(MAE_edge, group=='all', Duration=='first_5_mins'), #T = 10 min
                      val_name = 'rMSE_delta',
                      legTitle = '% Change\nin FC Error',
                      title = "% Change vs. No Censoring",
                      FD = ss,
                      zlim = c(-10, 10),
                      colFUN = seq3div))
    dev.off()
  }
  
  #quantify % of edges 
  summary <- MAE_edge %>%
    filter(group == "all", Duration == "first_5_mins") %>%
    group_by(FD) %>%
    summarise(
      above10 = 100 * sum(rMSE_delta > 10, na.rm = TRUE) / n(),
      below10 = 100 * sum(rMSE_delta < -10, na.rm = TRUE) / n(),
      above15 = 100 * sum(rMSE_delta > 15, na.rm = TRUE) / n(),
      below15 = 100 * sum(rMSE_delta < -15, na.rm = TRUE) / n(),
      above50 = 100 * sum(rMSE_delta > 50, na.rm = TRUE) / n(),
      below50 = 100 * sum(rMSE_delta < -50, na.rm = TRUE) / n()
    )
  
  write.table(summary, file = file.path(dir_github, "plots", "MAE", base, "dMAE_percent_edges.txt"), row.names = FALSE)
  

  ##########################################################################
  ### DURATION CHANGE

  load(file = file.path(dir_results, 'results','7_MAE',paste0('minDuration_edge_',base,'.RData'))) #MAE_minDur_edge

  MAE_minDur_edge$FD <- format_FD(MAE_minDur_edge$FD)
  MAE_minDur_edge$Duration_change_pct <- MAE_minDur_edge$Duration_change_pct*100

  for(ss in FD_levels_long[3:4]){
    fname_ss <- paste0('durationChange_', ss, '.png')
    png(file.path(dir_github, "plots", "MAE", base, fname_ss), width=1000, height=900, res=200)
    print(plot_FC_mat(df = filter(MAE_minDur_edge, group=='all'),
                val_name = 'Duration_change_pct',
                legTitle = '% Change\nin Duration',
                title = "% Change vs. FD5",
                FD = ss,
                zlim = c(-130, 130),
                colFUN = seq3div))
    dev.off()
  }
  
  MAE_minDur_edge %>%
    filter(!FD %in% c("None", "Lenient")) %>%
    group_by(FD, group) %>%
    summarize(min = min(Duration_change_pct, na.rm = TRUE),
              max = max(Duration_change_pct, na.rm = TRUE),
              over100 = mean(Duration_change_pct > 100, na.rm = TRUE))
  # P36: 
  # FD        group   min   max over100
  #   1 Stringent all   -77.1  500.  0.0783
  # 2 Stringent himo  -77.1  499.  0.102 
  # 3 Stringent lomo  -77.1  499.  0.0469
  # 4 Expanded  all   -77.1  500.  0.135 
  # 5 Expanded  himo  -77.1  500.  0.167 
  # 6 Expanded  lomo  -77.1  500.  0.0956
  # FIX:
  #   FD        group   min   max over100
  #   1 Stringent all   -77.1  491.  0.0441
  # 2 Stringent himo  -77.1  487.  0.0776
  # 3 Stringent lomo  -77.1  500.  0.0290
  # 4 Expanded  all   -77.1  491.  0.151 
  # 5 Expanded  himo  -77.1  499.  0.227 
  # 6 Expanded  lomo  -77.1  490.  0.0749
  

  ##########################################################################
  ### BASELINE NOISE

  ### BOXPLOTS (OVER EDGES) OF BASELINE NOISE (CHANGE VS NO SCRUBBING)

  load(file = file.path(dir_results, 'results','7_MAE',paste0('baselineSD_edge_',base,'.RData'))) #baselineSD_edge

  plot_baselineSD <- function(df, grp, colors, ylab){
    ggplot(filter(df, FD != "None", group==grp),
           aes(x = FD, y = baselineSDeff_mean_delta, fill = FD)) +
      geom_boxplot(outliers = FALSE) + #, color = 'black') +
      scale_fill_manual(values = colors) +
      #geom_point(size=2) +
      scale_y_continuous(labels = scales::percent) +#, limits = ylimit) +
      ggtitle('Baseline Noise vs. No Censoring') + ylab(ylab) +
      cowplot::theme_cowplot() +
      theme(legend.position = "none",
            axis.text.x=element_text(angle=45,hjust=1),
            axis.title.x = element_blank())
  }

  baselineSD_edge <- filter(baselineSD_edge, FD %in% FD_levels)
  baselineSD_edge$FD <- format_FD(baselineSD_edge$FD)

  pdf(file.path(dir_github, "plots", "MAE", base, "baselineSD_boxplot_eff_edge.pdf"), height=4, width=3)
  #below-average movers
  print(plot_baselineSD(df = baselineSD_edge, grp = 'Below-Average Motion', colors=myCols[-1], ylab="Change vs. No Censoring\n(Below-Average Motion Subjects)"))
  #above-average movers
  print(plot_baselineSD(df = baselineSD_edge, grp = 'Above-Average Motion', colors=myCols[-1], ylab="Change vs. No Censoring\n(Above-Average Motion Subjects)"))
  dev.off()

  ### MATRIX PLOTS OF BASELINE NOISE (CHANGE VS NO SCRUBBING)

  baselineSD_edge$baselineSDeff_mean_delta <- baselineSD_edge$baselineSDeff_mean_delta*100

  for(ss in FD_levels_long[2:4]){
    fname_ss <- paste0('baselineSD_', ss, '_himo.png')
    png(file.path(dir_github, "plots", "MAE", base, fname_ss), width=1000, height=900, res=200)
    print(plot_FC_mat(df = filter(baselineSD_edge, group=='Above-Average Motion'),
                      val_name = 'baselineSDeff_mean_delta',
                      legTitle = '% Change\nin Baseline Noise',
                      title = "% Change vs. No Censoring",
                      FD = ss,
                      zlim = c(-50, 50),
                      colFUN = seq3div))
    dev.off()
  }
}


