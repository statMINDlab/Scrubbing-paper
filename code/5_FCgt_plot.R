####################################################################
# Who: Joanne Hwang
# When: Spring 2026
# What: Create matrix visualizations for ground truth FC
# How: First, 5_error.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## -----------------------------------------------------------------------------
source("0_SharedCode.R")

source(file.path(dir_github, "code", "FC_vis_funs.R")) #tools and functions for FC matrix visualization

networklabels <- FALSE

plot_FCgt_mat <- function(df, legTitle, title, zlim = c(-15, 15), colFUN = seq3div){
  vals <- df
  
  if (is.data.frame(vals) && ncol(vals) == 1) {
    vals <- as.numeric(vals[[1]])
  }
  
  vals[vals < zlim[1]] <- zlim[1]
  vals[vals > zlim[2]] <- zlim[2]
  plt <- ggcorrplot2(
    values = vals,
    colFUN = colFUN,
    divColor="black",
    title = paste0(title),
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

for (base in c("P36","FIX_GSR")) {
  
  print(paste0('~~~~~~~~~~~~~~~~~~~~~~~~~~ ',base, ' ~~~~~~~~~~~~~~~~~~~~~~~~~~ '))
  
  FC_gt_edge <- readRDS(file = file.path(dir_github, 'results', '5_error', paste0('FCgt_over_subj_',base,'.rds')))
  
  png(file.path(dir_github, 'plots', paste0('FCgt_mat_',base,'.png')), width=1000, height=900, res=200)
  print(plot_FCgt_mat(df = FC_gt_edge,
                      legTitle = 'Ground Truth FC',
                      title = 'Ground Truth FC',
                      zlim = c(-1.0, 1.0),
                      colFUN = seq3div))
  dev.off()
}
