####################################################################
# Who: Amanda Mejia
# When: Fall 2024
# What: Plot QC-FC estimates and corrected p-values for each scrubbing level
# How: First, 7b_Analysis_QCFC_compute.R must be run
# Where: Can be run on HPC or locally, since it only requires access to Github files
####################################################################

## ---------------------------------------------------------------------------
source("0_SharedCode.R")

networklabels <- FALSE

#QCFC_naive_df --> cor(FD, FC)
#QCFC_adj_df --> partition QC-FC correlations into within- and between-subject associations using multiple linear regression
#QCFC_within_df --> compute within-subject associations in a simpler way as a sanity check

#loop over base denoising methods to make one big data frame
QCFC_naive_df0 <- QCFC_adj_df0 <- QCFC_within_df0 <- NULL
for (bb in 1:nB) {

  baseName <- baseNames[bb]

  load(file = file.path(dir_github, 'results', '7_QCFC', paste0('QCFC_df_',baseName,'.RData'))) #QCFC_naive_df, QCFC_adj_df, QCFC_within_df

  #remove FD4 and FD3
  QCFC_naive_df <- subset(QCFC_naive_df, !scrub %in% c('FD4', 'FD3'))
  QCFC_adj_df <- subset(QCFC_adj_df, !scrub %in% c('FD4', 'FD3'))
  QCFC_within_df <- subset(QCFC_within_df, !scrub %in% c('FD4', 'FD3'))

  QCFC_naive_df0 <- rbind(QCFC_naive_df0, QCFC_naive_df)
  QCFC_adj_df0 <- rbind(QCFC_adj_df0, QCFC_adj_df)
  QCFC_within_df0 <- rbind(QCFC_within_df0, QCFC_within_df)

} #end loop over base denoising methods

QCFC_naive_df <- QCFC_naive_df0
QCFC_adj_df <- QCFC_adj_df0
QCFC_within_df <- QCFC_within_df0
rm(QCFC_naive_df0, QCFC_adj_df0, QCFC_within_df0)

### PLOT TYPE 1: QC-FC in FC matrix form

source(file.path(dir_github, "code", "FC_vis_funs.R")) #tools and functions for FC matrix visualization

#this function loops over scrubbing levels (typically to make a multi-page PDF)
plot_fun <- function(df, val_name, baseName, p_adj = NULL, legTitle, title, FD_levels){
  #p_adj = method for p.adjust if applicable

  #loop over scrubbing (reordered from least to most)
  for (FD_ss in FD_levels) {

    #grab values from DF
    #FD_ss <- FD_levels[ss]
    vals <- df[,val_name]
    vals_bs <- vals[df$base==baseName & df$scrub==FD_ss]

    # if(!is.null(p_adj)){
    #   #for plotting p-values
    #   vals_bs <- p.adjust(vals_bs, method = p_adj)
    #   vals_bs <- 1*(vals_bs < 0.01)
    #   lims <- c(0, 1)
    #   colFUN = seq2
    # } else {
      #for plotting QC-FC correlation estimates
    lims <- c(-0.3, 0.3)
    colFUN = seq3div2 #updated --> (turquoise, white, orange)
    # }

    vals_bs[vals_bs < lims[1]] <- lims[1]
    vals_bs[vals_bs > lims[2]] <- lims[2]
    plt <- ggcorrplot2(
      values = vals_bs,
      colFUN = colFUN,
      divColor="black",
      title = paste0(title, ", ", FD_ss),
      legTitle = legTitle,
      lim=lims
    ) +
      theme_classic() +
      theme(
        legend.position = "bottom",
        axis.title.x = element_text(size = 36, margin = margin(t = 10)),
        axis.title.y = element_text(size = 36, margin = margin(t = 10)),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) +
      labs(x = "Node", y = "Node")

    print(plt)
  }
}

### Visualize matrices of QC-FC measures

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#FD_reo <- FD_levels2[-4]
nFD <- length(FD_levels)

#loop over denoising
for (bb in 1:nB) {

  print(baseName <- baseNames[bb])

  ### 1) Plot QC-FC estimates for each scrubbing level
  ### 2) Plot corrected p-values for each scrubbing level


  #loop over Naive, Within-Subject and Between-Subject QCFC
  for(meas in c('naive','within','between')){

    df <- if(meas == 'naive') QCFC_naive_df else QCFC_adj_df
    val_name <- if(meas == 'naive') "QCFC" else paste0("QCFC_", meas)
    title <- if(meas == 'naive') meas else paste0(meas,'-Subject')
    title <- paste0(firstup(title), ' QC-FC')

    #PDFs (all scrubbing levels)

    fname <- paste0('QCFC_',firstup(meas),'.pdf')
    pdf(file.path(dir_github, "plots", "QCFC", baseName, fname))
    suppressMessages(plot_fun(df = df, val_name = val_name, baseName = baseName,
                              legTitle = 'QC-FC', title = title, FD_levels = FD_levels))
    dev.off()

    #PNGs (one per scrubbing level)

    #nS <- length(FD_reo)
    for(ss in 1:nFD){

      FD_ss <- FD_levels[ss]

      df <- if(meas == 'naive') QCFC_naive_df else QCFC_adj_df
      val_name <- if(meas == 'naive') "QCFC" else paste0("QCFC_", meas)
      title <- if(meas == 'naive') meas else paste0(meas,'-Subject')
      title <- paste0(firstup(title), ' QC-FC')

      fname <- paste0('QCFC_',firstup(meas),'_',FD_ss,'.png')
      png(file.path(dir_github, "plots", "QCFC", baseName, fname), width=600, height=600)
      suppressMessages(plot_fun(df = df, val_name = val_name, baseName = baseName,
                                legTitle = 'QC-FC', title = title, FD_levels = FD_ss))
      dev.off()


    } #end loop over scrubbing levels

  } #end loop over QC-FC measures

  #PNGs of alternative Within-Subject QCFC

  for(grp in c('all','hivar')){

    for(ss in 1:nFD){

      FD_ss <- FD_levels[ss]

      fname <- paste0('QCFC_Within2_',grp,'_',FD_ss,'.png')
      png(file.path(dir_github, "plots", "QCFC", baseName, fname), width=600, height=600)
      suppressMessages(plot_fun(df = QCFC_within_df, val_name = grp, baseName = baseName,
                                legTitle = 'QC-FC', title = paste0("Within-Subject QC-FC (", grp,")"), FD_levels = FD_ss))
      dev.off()
    }
  }

} # end loop over baseline denoising methods



# HERE -- for each edge, average QC-FC measures among distal nodes for each node and make images (focus on between-subject)


# ### PLOT: QC-FC averaged across connections, plotted on the brain
#
# plot_fun <- function(xii, fname_start, fname_end, zlim = c(-0.35, 0.35), titles = unique(Labs_reo$network2)[1:8], fnames = Labs_networks3){
#   #pal <- c("white","gray90", "pink1", "red","red4")
#   pal <- c("darkblue","dodgerblue","skyblue","gray90","white","gray90", "pink1", "red","red4")
#   plot(xii, zlim = zlim,
#        colors=pal,
#        idx = 1:ncol(xii),
#        title = titles,
#        fname = paste0(fname_start, fnames, fname_end))
# }
#
# #everything except DMN, FP (too close), and L (inconsistent info on what it is)
# inds_nonDMN <- c(1,4,5:8)
#
# for (bb in 3:nB) {
#   print(baseName <- baseNames[bb])
#
#   for (ss in 1:nScrub) {
#
#     print(FD_ss <- scrubNames2[ss]) #partial file name
#     fname_end <- paste0("_", FD_ss)
#     fname_start <- file.path('plots','QCFC',baseName, 'QCFC_')
#
#     #Naive QCFC
#     vals_bs <- QCFC_naive_df$QCFC[QCFC_naive_df$base==baseName & QCFC_naive_df$scrub==FD_ss]
#     xii_bs <- ggImgplot2(vals_bs, zlim = c(-0.35, 0.35))
#     plot_fun(xii_bs, fname_start = paste0(fname_start, "Naive_"), fname_end = fname_end)
#
#     #Within-Subject QCFC
#     vals_bs <- QCFC_adj_df$QCFC_within[QCFC_adj_df$base==baseName & QCFC_adj_df$scrub==FD_ss]
#     xii_bs <- ggImgplot2(vals_bs, zlim = c(-0.35, 0.35))
#     plot_fun(xii_bs, fname_start = paste0(fname_start, "Within_"), fname_end = fname_end)
#     xii_bs2 <- newdata_xifti(xii_bs, rowMeans(as.matrix(xii_bs)[,inds_nonDMN]))
#     plot_fun(xii_bs2, titles = "nonDMN", fnames="nonDMN", fname_start = paste0(fname_start, "Within_"), fname_end = fname_end, zlim = c(-0.35, 0.35))
#
#     #Between-Subject QCFC
#     vals_bs <- QCFC_adj_df$QCFC_between[QCFC_adj_df$base==baseName & QCFC_adj_df$scrub==FD_ss]
#     xii_bs <- ggImgplot2(vals_bs, zlim = c(-0.35, 0.35))
#     plot_fun(xii_bs, fname_start = paste0(fname_start, "Between_"), fname_end = fname_end)
#     xii_bs2 <- newdata_xifti(xii_bs, rowMeans(as.matrix(xii_bs)[,inds_nonDMN]))
#     plot_fun(xii_bs2, titles = "nonDMN", fnames="nonDMN", fname_start = paste0(fname_start, "Between_"), fname_end = fname_end, zlim = c(-0.35, 0.35))
#
#   }
# }
#
#