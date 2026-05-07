####################################################################
# Who: Joanne Hwang and Amanda Mejia
# When: Fall 2024
# What: Plot QC-FC measures as a function of Euclidean distance
# How: First, 7b_Analysis_QCFC_compute.R must be run to compute QCFC and
# construct data frames
# Where: Can be run on HPC or locally, since it only requires access to Github files
####################################################################

## ---------------------------------------------------------------------------
source("0_SharedCode.R")

#load additional packages
library(ggthemes)
library(ggpointdensity)
library(viridis)
library(reshape2)
library(dplyr)

#create vector with reordering of parcels
reordering <- read.csv(file.path(dir_github, "data", "reordering.csv"))
reordering_v <- pull(reordering, A) #A indicates the first column of the parcels (of A/S/C)
print(reordering_v) #resulting vector with roi reordering

#read in the csv file
Schaefer_Yeo17_xyz <- read.csv(file.path(dir_github, "data", "Parc", "Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv"))

#reorder the rows of the table
Schaefer_Kong17_xyz <- Schaefer_Yeo17_xyz[order(reordering_v),] #Schaefer_Yeo17 parcels reordered to be in Schaefer_Kong17 order

#compute euclidean distance
euclidean_distance <- as.matrix(dist(Schaefer_Kong17_xyz[,3:5], method = "euclidean"))

#visualize distance matrix
pdf(file.path(dir_github, "plots", "DistanceMatrix.pdf"), width = 6, height = 6)
plot_FC(euclidean_distance, zlim = c(-200, 200))
dev.off()

#vectorize upper triangle of distance matrix
ut_vector <- euclidean_distance[upper.tri(euclidean_distance)] #has 79800 elements

## Load QCFC dataframes ---------------------------------------------------------------------------

#vector of inclusion in 400x400 (these are the edges for which we can compute)
mat <- matrix(0, 419, 419) #create 419x419 matrix containing all 0 values
mat[1:400,1:400] <- 1 #columns 1-400 and rows 1-400 contain value of 1
mat_UT <- mat[upper.tri(mat)] #vectorize upper triangle (0.5*(400*399) = 319200 elements)
nEdge400 <- 400*399/2

## ---------------------------------------------------------------------------

#PLOT QCFC BY DISTANCE

#loop over base denoising methods
nFD <- length(FD_levels) #FD_levels <- c("None","FD5","FD2","FD2+")
#nFD <- length(FD_levels2) #including FD3 and FD4

#plot QCFC measures by euclidean distance using QCFC_all_df (QCFC_adj_df and QCFC_naive_df)
for(bb in baseNames){
  # if (bb=="P32") { next } # DAMON TEMP

  print(bb)

  # load in QCFC data frames

  #QCFC_naive_df, QCFC_adj_df, QCFC_within_df
  load(file = file.path(dir_github, 'results', '7_QCFC', paste0('QCFC_df_',bb,'.RData')))

  #remove FD4 and FD3
  QCFC_naive_df <- subset(QCFC_naive_df, !scrub %in% c('FD4', 'FD3'))
  QCFC_adj_df <- subset(QCFC_adj_df, !scrub %in% c('FD4', 'FD3'))
  QCFC_within_df <- subset(QCFC_within_df, !scrub %in% c('FD4', 'FD3'))

  #add node-to-node distances to QCFC data frames
  QCFC_adj_df$in400 <- QCFC_naive_df$in400 <- QCFC_within_df$in400 <- mat_UT #create new column called 'in400' with 1 values from mat_UT
  QCFC_adj_df <- subset(QCFC_adj_df, in400==1) #only include edges in 400x400 submatrix, has 1276800 elements (319200*4)
  QCFC_naive_df <- subset(QCFC_naive_df, in400==1)
  QCFC_within_df <- subset(QCFC_within_df, in400==1)
  QCFC_adj_df$dist <- QCFC_naive_df$dist <- QCFC_within_df$dist <- ut_vector #add in distance between pairs of 400 cortical nodes
  #QCFC_adj_df$scrub <- QCFC_naive_df$scrub <- QCFC_within_df$scrub <- factor(QCFC_adj_df$scrub, levels = FD_levels) #categorize by level of scrubbing
  QCFC_adj_df$scrub <- QCFC_naive_df$scrub <- QCFC_within_df$scrub <- factor(QCFC_adj_df$scrub, levels = FD_levels2) #categorize by level of scrubbing


  ### PLOT TYPE 1a: Plot naive and adjusted QCFC by euclidean distance

  #add naive QCFC to df of adjusted QCFC for comparison
  QCFC_adj_df$QCFC_naive <- QCFC_naive_df$QCFC

  #melt into long format
  QCFC_adj_df$edge <- rep(1:nEdge400, nFD) #df grouped by scrubbing methods
  QCFC_long <- melt(QCFC_adj_df[,c('QCFC_within','QCFC_between','QCFC_naive','scrub','dist','edge')],
                    id.vars = c('edge','scrub','dist'))

  QCFC_long$variable <- factor(QCFC_long$variable,
                              levels = c('QCFC_naive', 'QCFC_between', 'QCFC_within'),
                              labels = c('Naive QC-FC','Between-Participant QC-FC','Within-Participant QC-FC'))

  QCFC_long$scrub <- factor(QCFC_long$scrub,
                            levels = FD_levels,
                            labels = FD_levels_long)

  QCFC_adj_dist_cor <- summarize(group_by(QCFC_long, scrub, variable, .groups='keep'),
                                    mycor = cor(value, dist),
                                    count = n())
  QCFC_adj_dist_cor$mycor2 <- paste0('r = ', format(round(QCFC_adj_dist_cor$mycor, 3), nsmall=3))

  QCFC_types <- levels(QCFC_long$variable)

  for (v in QCFC_types) {

    yscale <- if (bb == "P36") { c(-0.1,0,0.1) } else { seq(-0.1,0.3,by=0.1) }
    yaxis <- list("Within-Participant QC-FC" = list(ytitle = "Within-Subject QC-FC"), # Repeated Measures QC-FC\n (Motion Artifact/Slate)
                  "Between-Participant QC-FC" = list(ytitle = "Repeated Measures QC-FC\n (Motion Trait)"),
                  "Naive QC-FC" = list(ytitle = "Standard QC-FC"))

    df_v <- subset(QCFC_long, variable == v)
    cor_v <- subset(QCFC_adj_dist_cor, variable == v)
    fname <- paste0(gsub("[^A-Za-z0-9]+", "_", v), "_Adjusted_Distance.pdf")

    custom_labels <- as_labeller(c(
      "None" = "No Scrubbing",
      "Lenient" = "Lenient",
      "Stringent" = "Stringent",
      "Expanded" = "Expanded"
    ))

    pdf(file.path(dir_github, "plots", "QCFC", bb, fname), width=7, height=4) #QCFC_distance
    print(ggplot(df_v, aes(x=dist, y=value, group=scrub)) +
            stat_density_2d(geom='polygon', alpha = 0.2, fill='royalblue') +
            geom_hline(yintercept=0, col='black') +
            geom_smooth(method="lm", col='red') +
            geom_text(data = cor_v, aes(label = mycor2, x=Inf, y=Inf, hjust=1.15, vjust=1.5), col='red', size = 4) +
            scale_y_continuous(breaks = yscale) +
            #theme_few() +
            cowplot::theme_cowplot() +
            facet_grid(col = vars(scrub), scales='free', labeller=custom_labels) +
            labs(title = v, x = 'Node-Node Distance', y = yaxis[[v]]$ytitle) +
            theme(legend.position='none',
                  strip.background = element_blank(),
                  strip.text = element_text(face = 'bold')))
    dev.off()
  }

  ### PLOT TYPE 1b: Plot trend lines for naive and adjusted QCFC by euclidean distance

  #add a fake variable for faceting
  #QCFC_long$blank <- ""

  myCols <- c("black", "#5585bd", "#b3ce91", "#b63545") #for None / FD5 / FD2 / FD2+
  yscale <- if (bb == "P36") { c(0,0.02,0.04) } else { seq(0,0.40,by=0.1) }
  ycoord <- if (bb == "P36") { c(-0.01,0.045) } else { c(-0.04,0.40) }

  for (v in QCFC_types) {

    yaxis <- list("Within-Participant QC-FC" = list(ytitle = "Repeated Measures QC-FC\n (Motion Artifact/Slate)"),
                  "Between-Participant QC-FC" = list(ytitle = "Repeated Measures QC-FC\n (Motion Trait)"),
                  "Naive QC-FC" = list(ytitle = "Standard QC-FC"))

    df_v <- subset(QCFC_long, variable == v)
    fname <- paste0(gsub("[^A-Za-z0-9]+", "_", v), "_Adjusted_Distance_Trend.pdf")
    fname2 <- paste0(gsub("[^A-Za-z0-9]+", "_", v), "_Adjusted_Distance_Trend_Smooth.pdf")

    pdf(file.path(dir_github, "plots", "QCFC", bb, fname), width=3.5, height=4) #QCFC_distance_trend
    print(ggplot(df_v, aes(x=dist, y=value, color=scrub, group=scrub)) +
            geom_hline(yintercept=0, col='black') +
            geom_smooth(method='lm', se=FALSE, linewidth = 0.7) +
            scale_color_manual(values = myCols) +
            scale_y_continuous(breaks = yscale) +
            coord_cartesian(ylim = ycoord) +
            #theme_few() +
            cowplot::theme_cowplot() +
            #facet_grid(rows=vars(blank), col= vars(scrub), scales='free') +
            labs(title = paste0(v, "\n"), x='Node-Node Distance', y=yaxis[[v]]$ytitle) +
            theme(legend.position=c(0.97,0.98),
                  legend.justification=c("right","top"),
                  legend.direction = "vertical",
                  legend.title=element_blank(),
                  legend.text=element_text(size=8.5),
                  legend.key.size = unit(0.5, "cm"),
                  legend.spacing.y = unit(0.1, "cm")))
    dev.off()

    # pdf(file.path(dir_github, "plots", "QCFC", bb, fname2), width=4, height=3.5) #QCFC_distance_trend_smooth
    # print(ggplot(df_v, aes(x=dist, y=value, color=scrub, group=scrub)) +
    #         geom_hline(yintercept=0, col='black') +
    #         geom_smooth(method='loess', se=FALSE, linewidth = 0.7) +
    #         scale_color_manual(values = myCols) +
    #         scale_y_continuous(breaks = yscale) +
    #         #coord_cartesian(ylim = ycoord) +
    #         theme_few() +
    #         #facet_grid(rows=vars(blank), col= vars(scrub), scales='free') +
    #         labs(title = v, x='Node-Node Distance', y='Motion-FC Correlation') +
    #         theme(legend.position='bottom',
    #               legend.title=element_blank(),
    #               legend.text=element_text(size=8.5),
    #               legend.spacing.x=unit(0.05, "cm")))
    # dev.off()
  }

  ### PLOT TYPE 2: Plot simple within-subject QCFC by euclidean distance

  QCFC_within_df$edge <- rep(1:nEdge400, nFD)

  pdf(file.path(dir_github, "plots", "QCFC", bb, "QCFC_Within_Distance.pdf"), width=8, height=3.5) #QCFC_distance
  #loop over subsets of subjects
  for(s in c('all','hivar')){ #,'hivar_hiFD','hivar_loFD')

    QCFC_within_dfs <- QCFC_within_df[,c(s,'scrub','dist','edge')]
    names(QCFC_within_dfs)[1] <- 'QCFC'

    QCFC_within_dfs$scrub <- factor(QCFC_within_dfs$scrub,
                              levels = FD_levels,
                              labels = FD_levels_long)

    QCFC_within_dist_cor <- summarize(group_by(QCFC_within_dfs, scrub),
                               mycor = cor(QCFC, dist),
                               count = n())
    QCFC_within_dist_cor$mycor2 <- paste0('r = ', format(round(QCFC_within_dist_cor$mycor, 3), nsmall=3))

    #remove FD4
    # QCFC_within_dfs <- subset(QCFC_within_dfs, scrub != 'FD4')
    # QCFC_within_dist_cor <- subset(QCFC_within_dist_cor, scrub != 'FD4')

    print(ggplot(QCFC_within_dfs, aes(x=dist, y=QCFC)) +
            stat_density_2d(geom='polygon', alpha = 0.2, fill='royalblue') +
            geom_hline(yintercept=0, col='black') +
            geom_smooth(method="lm", col='red') +
            geom_text(data = QCFC_within_dist_cor, aes(label = mycor2, x=Inf, y=Inf, hjust=1.15, vjust=1.5), col='red') +
            theme_few() + facet_grid(. ~ scrub, scales='free') +
            labs(title = paste0("QCFC ", s), x = 'Node-Node Distance', y = 'Motion-FC Correlation') +
            theme(legend.position='none')) #plot dist by QCFC grouped by scrubbing level for each cleaning method
  }
  dev.off()

} #end loop over baseline denoising

