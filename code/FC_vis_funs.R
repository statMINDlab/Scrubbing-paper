####################################################################
# Who: Amanda Mejia?
# When: Fall 2024?
# What: Parcellation information for FC matrix visualization
# How: First, 1_ParcMatLabels.R must be run
# Where: Can be run on RED or locally, since it only requires access to Github files
####################################################################

## ---------------------------------------------------------------------------
#Parcellation
Labs <- readRDS(file.path(dir_github, "data/Parc/ParcLabels.rds"))$Labs
Labs$label2 <- gsub("_.*", "", gsub("-L|-R", "", gsub("RH_|LH_", "", Labs$label)))
Labs$label2 <- factor(Labs$label2, levels=unique(Labs$label2))
Labs_reo <- Labs[order(Labs$idx2),]
Labs_networks <- levels(unique(Labs_reo$label2))
Labs_xii <- read_cifti(file.path(dir_github, "data/Parc", "Schaefer2018_400Parcels_Kong2022_17Networks_order.dlabel.nii"))

#Parcel coordinates
xyz <- read.csv(file.path(dir_github, "data/Parc/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv"))
xyz$ROI.Name <- gsub('17Networks_','',xyz$ROI.Name)

# Variables
mbar_cols <- ifelse(is.na(Labs$group_color), Labs$network_color, Labs$group_color)[order(Labs$idx2)]
cor_mat_ylabs <- ifelse(
  c(Labs$network_first[rev(order(Labs$idx2))][seq(3, parc_res2)], FALSE, FALSE),
  Labs$network3[rev(order(Labs$idx2))], ""
)
gg_pdv <- which(Labs$network_first[order(Labs$idx2)])
gg_pdv <- gg_pdv[gg_pdv != 1]
divwidth <- .3
gg_pdv <- data.frame(xmin=gg_pdv-.5-(divwidth/2), xmax=gg_pdv-.5+(divwidth/2))
rm(divwidth)
gg_pdv$Var1 <- Var2 <- 0
pdv <- c(.5, gg_pdv$xmin + .15, parc_res2+.5) + .5 # "parc_div"


# Converted a vectorized lower triangular correlation matrix back to its full form.
cor_mat <- function(x_diag, diag_val=NA, names=NULL, newOrder=NULL, lowerOnly=FALSE) {
  d <- 1/2 + sqrt(1/4 + 2*length(x_diag))
  mat <- diag(d) * diag_val
  mat[upper.tri(mat)] <- x_diag
  mat <- t(mat)
  mat[upper.tri(mat)] <- x_diag
  if (!is.null(names)) {
    stopifnot(length(names)==d); rownames(mat) <- colnames(mat) <- names
  }
  if (!is.null(newOrder)) {
    stopifnot(all(sort(newOrder) == seq(d))); mat <- mat[newOrder,rev(newOrder)]
  }
  if (lowerOnly) { mat[seq(d), rev(seq(d))][lower.tri(mat)] <- NA }
  mat
}



ggcorrplot2 <- function(values, colFUN, title, legTitle, divColor="black", lim=NULL, diagVal=0, legend.position='bottom') {
  # Get cor mat
  mat <- values
  mat <- fMRItools:::cor_mat(mat, names=Labs$label, newOrder=order(Labs$idx2))

  # # Take network-network means in upper tri
  # mat2 <- mat
  # for (jj in seq(nrow(gg_pdv)+1)) {
  #   for (kk in seq(nrow(gg_pdv)+1)) {
  #     mat2[seq(pdv[jj], pdv[jj+1]-1), parc_res2+1 - seq(pdv[kk], pdv[kk+1]-1)] <- mean(
  #       mat2[seq(pdv[jj], pdv[jj+1]-1), parc_res2+1 - seq(pdv[kk], pdv[kk+1]-1)], na.rm=TRUE
  #     )
  #   }
  # }
  mat[is.na(mat)] <- diagVal
  # mat[seq(nrow(mat)), rev(seq(ncol(mat)))][lower.tri(mat)] <- mat2[seq(nrow(mat)), rev(seq(ncol(mat)))][lower.tri(mat)]

  # Limits
  if (!is.null(lim)) {
    if (length(lim)==1) { lim <- c(-lim, lim) }
    mat[] <- pmax(lim[1], pmin(lim[2], mat))
  } else {
    lim <- max(abs(mat[]))
    lim <- c(-lim, lim)
  }

  p <- ggcorrplot::ggcorrplot(mat, outline.color = "#00000000", title=title, digits=12) +
    scale_y_discrete(labels=cor_mat_ylabs) +
    #coord_equal(xlim=c(-4-10, parc_res2+.65), ylim=c(.5, parc_res2+5+10), expand=FALSE, clip = "off") +
    #coord_equal(xlim=c(0.5, parc_res2+.05), ylim=c(0.5, parc_res2+0.5), expand=FALSE, clip = "off") +
    colFUN(
      c(lim[1], lim[2]), guide=guide_colorbar(ticks.colour = divColor, ticks.linewidth = 1),
      labels = function(x){gsub("0.", ".", x, fixed=TRUE)}, na.value="yellow"
    ) +
    labs(fill=legTitle) +
    theme(
      plot.margin = unit(c(5,5,5,5), "pt"),
      panel.grid.major = element_blank(),
      axis.text.y = element_text(margin=margin(r=10)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = legend.position,
    )

  if (networklabels == TRUE) {
    p <- p + coord_equal(xlim=c(-4-10, parc_res2+.65),
                         ylim=c(.5, parc_res2+5+10),
                         expand=FALSE, clip = "off")
  } else {
    p <- p + coord_equal(xlim=c(0.5, parc_res2+.05),
                         ylim=c(0.5, parc_res2+0.5),
                         expand=FALSE, clip = "off")
    }

  if (networklabels == TRUE) {
    p <- p +

      #brain region label color bar
      annotate("rect", xmin=-4-10, xmax=.5, ymin=parc_res2+.5-seq(nrow(mat)), ymax=parc_res2+.5-seq(nrow(mat))+1,
               fill=mbar_cols[seq(nrow(mat))]) +
      annotate("rect", ymin=parc_res2+.5, ymax=parc_res2+5+10, xmin=seq(nrow(mat))-.5, xmax=seq(nrow(mat))+.5,
               fill=mbar_cols[seq(nrow(mat))]) +
      # annotate("rect", xmin=-4-2, xmax=.5, ymin=parc_res2+.5-seq_len(nrow(mat)), ymax=parc_res2+.5-seq_len(nrow(mat))+1,
      #          fill=mbar_cols[seq_len(nrow(mat))], color = NA) +
      # annotate("rect", ymin=parc_res2+.5, ymax=parc_res2+1.5, xmin=seq_len(nrow(mat))-.5, xmax=seq(nrow(mat))+.5,
      #          fill=mbar_cols[seq_len(nrow(mat))], color = NA) +

      #dividers
      geom_rect(aes(xmin=xmin, xmax=xmax, ymin=0, ymax=parc_res2+5+10), data=gg_pdv, fill=divColor) +
      geom_rect(aes(ymin=parc_res2+1-xmin, ymax=parc_res2+1-xmax, xmin=-4-10, xmax=parc_res2+.5), data=gg_pdv, fill=divColor) +

      # #four edges
      annotate("rect", xmin=.35, xmax=.65, ymin=0, ymax=parc_res2+5+10, fill=divColor) +
      annotate("rect", xmin=parc_res2+.35, xmax=parc_res2+.65, ymin=0, ymax=parc_res2+5+10, fill=divColor) +
      annotate("rect", ymin=.35, ymax=.65, xmin=-4-10, xmax=parc_res2+.85, fill=divColor) +
      annotate("rect", ymin=parc_res2+.35, ymax=parc_res2+.65, xmin=-4-10, xmax=parc_res2+.85, fill=divColor) +

      #separate brain region label color bar
      annotate("rect", ymin=0, ymax=parc_res2+.65, xmin=-.5-3, xmax=.35, fill="white") +
      annotate("rect", ymin=parc_res2+.65, ymax=parc_res2+1.5+3, xmin=-.5, xmax=parc_res2+.65, fill="white")
  } else {
    #add border around all edges
    p <- p + annotate("rect", xmin = 0.5, xmax = parc_res2 + 0.5,
             ymin = 0.5, ymax = parc_res2 + 0.5,
             color = "black", fill = NA, linewidth = 1)
  }
  return(p)
}

color_pals <- list(
  Beach2 = data.frame(
    c = rev(c(
      "#2c308c", "#4f52a3", "#6985f2", "#c7f1f9", "#e0fdec",
      "#fffff5", "#fefbcb", "#f7ea74", "#ee8718", "#d03a34", "#a81e22"
    )),
    v = 1 - c(0, .1, .225, .35, .425, .5, .575, .65, .775, .9, 1)
  )
)

seq1 <- function(limits=c(0,.1), ...){
  viridis::scale_fill_viridis(option="inferno", limits=limits, ...)
}

seq1_reverse <- function(limits=c(0,.1), ...){
  viridis::scale_fill_viridis(option="inferno", limits=limits, direction = -1, ...)
}

seq2 <- function(limits=c(0,.1), ...){
  scale_fill_gradient(low = 'white', high='purple', limits=limits, ...)
}

seq3 <- function(limits=c(0,.1), ...){
  scale_fill_gradientn(colors = c("white","gray90", "pink1", "red","red4"), limits=limits, ...)
}

seq3div <- function(limits=c(0,.1), ...){
  scale_fill_gradientn(colors = c("darkblue","dodgerblue","skyblue","gray90","white","gray90", "pink1", "red","red4"), limits=limits, ...)
}

seq3div2 <- function(limits=c(0,.1), ...){
  scale_fill_gradientn(colors = c("lightseagreen","turquoise","cadetblue2","gray90","white","gray90", "darksalmon", "orange","darkorange"), limits=limits, ...)
}

seq4 <- function(limits=c(0,.1), ...){
  viridis::scale_fill_viridis(option="viridis", limits=limits, ...)
}

seq5 <- function(limits=c(0,.1), ...){
  viridis::scale_fill_viridis(option="plasma", limits=limits, ...)
}

seq5_reverse <- function(limits=c(0,.1), ...){
  viridis::scale_fill_viridis(option="plasma", limits=limits, direction = -1, ...)
}

seq6div <- function(limits=c(0,.1), ...){
  vals <- brewer_pal(palette = 'PiYG', direction = -1)(9)
  scale_fill_gradientn(colors = vals, limits=limits, ...)
}

seq7 <- function(limits=c(0,.1), ...){
  viridis::scale_fill_viridis(option="rocket", limits=limits, direction = -1, ...)
}


#average across corr mat rows or columns to get an image

Labs_networks3 <- unique(Labs_reo$network3)
Labs_networks3 <- Labs_networks3[-length(Labs_networks3)] #remove SC for now

Labs2_xii <- Labs_xii
rownames(Labs2_xii$meta$cifti$labels[[1]]) <- gsub("17networks_","", rownames(Labs2_xii$meta$cifti$labels[[1]]))
Labs2_names <- rownames(Labs2_xii$meta$cifti$labels[[1]])[-1]
Labs2_mat <- as.matrix(Labs2_xii)

ggImgplot2 <- function(values, fname_start, fname_end, zlim = c(0.1, 0.35), highlightDMN = FALSE){
  # Get cor mat (419x419)
  mat <- values
  mat <- fMRItools:::cor_mat(mat, names=Labs$label, newOrder=order(Labs$idx2))
  #all.equal(rownames(mat)[1:400], Labs_reo$label[1:400]) #TRUE

  #loop over networks
  dat <- matrix(0*Labs2_mat, nrow=nrow(Labs2_mat), ncol=length(Labs_networks3))
  dat[Labs2_mat==0,] <- NA #medial wall locations

  for(n in 1:length(Labs_networks3)){
    net <- Labs_networks3[n]
    rows_n <- which(Labs_reo$network3 == net)
    mat_n <- mat[rows_n,]
    Img_n <- rev(colMeans(mat_n, na.rm=TRUE))
    #all.equal(names(Img_n), Labs_reo$label) #TRUE
    Img_n <- Img_n[1:400] #exclude SC for now
    if(highlightDMN){
      rowsDMN <- grepl("Default", names(Img_n))
      Img_n[!rowsDMN] <- (-1)*Img_n[!rowsDMN]
    }

    #replace xii data with values in Img_n
    for(k in 1:400){
      lab_k <- Labs2_names[k]
      inds_k <- which(Labs2_mat==k)
      val_k <- Img_n[names(Img_n) == lab_k]
      dat[inds_k,n] <- val_k
    }
  }

  if(highlightDMN){
    zlim <- c(-zlim[2], zlim[2])
    pal <- c("royalblue4","royalblue","lightskyblue1","gray90","white","gray90", "pink1", "red","red4")
  } else {
    pal <- c("white","gray90", "pink1", "red","red4")
  }

  newdata_xifti(convert_to_dscalar(Labs2_xii), dat)

}

