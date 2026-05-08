####################################################################
# Who: Amanda Mejia
# When: Fall 2024 to Spring 2025
# What: Obtain BWAS values and plot estimated vs. true correlation, attenuation, and mean bias
# How: First, run 8_BWAS_compute.R
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

source("0_SharedCode.R")
library(dplyr) #version 1.1.4
library(tidyr) #version 1.3.1
library(viridis) #version 0.6.5

### Read in df_rho and visualize true rho values
source('code/FC_vis_funs.R')
lims <- c(-0.15, 0.15)
lims_ICC <- c(0.5, 1)

networklabels <- FALSE

#####################################################################

# Read and visualize ICC of behavioral variables

ICC_demo <- readRDS(file = file.path(dir_github, 'results', '8_BWAS', 'ICC_demo.rds'))
ICC_demo <- sort(ICC_demo, decreasing = TRUE)
ICC_y <- ICC_demo[1] #maximal ICC value for any demographic variable (CogTotalComp_Unadj)
ICC_demo_df <- data.frame(variable = names(ICC_demo), ICC = ICC_demo)
ICC_demo_df$variable <- factor(ICC_demo_df$variable, levels = rev(ICC_demo_df$variable))
ICC_demo_df$ICC_group <- ifelse(ICC_demo_df$ICC > 0.9, "excellent",
                                ifelse(ICC_demo_df$ICC > 0.75, "good",
                                       ifelse(ICC_demo_df$ICC > 0.5, "moderate", "poor")))


pdf(file.path(dir_plots, "BWAS", "ICC_demo.pdf"), height=7, width=7)
ggplot(ICC_demo_df, aes(x = variable, y = ICC, fill=ICC_group)) +
  geom_hline(yintercept = c(0.5, 0.75, 0.9), lty=2, lwd=0.3, col='gray') +
  geom_bar(stat='identity') + coord_flip() +
  ylim(0, 1) + scale_fill_viridis_d(option = 'F') + theme_few() +
  theme(axis.text.y = element_text(size=8),
        axis.title.y= element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.87, 0.1),
        legend.text.position = 'left',
        legend.title = element_blank())
dev.off()

#####################################################################

### Visualize BWAS variance vs. N for low-reliability and high-reliability FC

N <- seq(500, 1500, 10)
rho <- 0.1
load(file = file.path(dir_github, 'results', '8_BWAS', paste0('BWAS_P36_ICC_FC.RData')))
(ICC_FC_quants <- quantile(ICC_X, probs = c(0.05, 0.5, 0.95), na.rm=TRUE))

var_rho_df <- expand.grid(N = N, ICC_FC = c(ICC_FC_quants, 1), ICC_y = ICC_y, var = NA)
var_rho_df$var <- (1/var_rho_df$N) * (1 - rho^2 * var_rho_df$ICC_FC * var_rho_df$ICC_y)^2
var_rho_df$ICC_FC_grp <- factor(var_rho_df$ICC_FC, levels = rev(c(ICC_FC_quants, 1)),
                                labels = rev(c('low (Q05)', 'medium (med)', 'high (Q95)', 'perfect')))

pdf(file.path(dir_github, 'plots', 'BWAS', 'BWAS_var_math.pdf'), width=6, height=4)
cols <- c('black','royalblue','darkorchid','deeppink')
hline_val <- (var_rho_df$var[var_rho_df$N==1500 & var_rho_df$ICC_FC_grp=='low (Q05)']) #BWAS variance for N=1500 for low ICC_FC
vline_val <- 1/hline_val * (1 - rho^2 * ICC_FC_quants[3] * ICC_y)^2 #theoretical value of N to achieve same BWAS variance for high ICC_FC
ggplot(var_rho_df, aes(x = N, y = sqrt(var), group= ICC_FC_grp, color=ICC_FC_grp)) +
  geom_line() +
  geom_vline(xintercept = 1500, lty=2, col='deeppink', alpha = 0.5) +
  geom_hline(yintercept = sqrt(hline_val), lty=2, alpha = 0.5) +
  geom_vline(xintercept = vline_val, lty=2, col='royalblue', alpha = 0.5) +
  ylab('BWAS Standard Deviation') + xlab('Sample Size') +
  scale_color_manual('ICC of FC', values = (cols)) +
  theme_few()
dev.off()

################################################################

### Read in BWAS results and ICC of FC

df_rho <- NULL
BWAS_propBias_math_mean <- list()
for (bb in 2:nB) {

  (baseName <- baseNames[bb])

  df_rho_bb <- NULL
  #read in df_rho
  for (ss in c(1,2,5,6)) {
    pfname_ss <- scrubNames[ss] #partial file name
    scrub_ss <- FD_levels2[ss]
    df_rho_bs <- readRDS(file=file.path(dir_github, 'results', '8_BWAS', paste0('BWAS_',baseName,'_',FD_levels2[ss],'.rds')))
    df_rho_bb <- rbind(df_rho_bb, df_rho_bs)
  } #end loop over scrubbing levels

  #visualize true value of rho
  pdf(file.path(dir_plots, "BWAS", paste0("rho_true_", baseName,".pdf")), height=6, width=6)
  rho_true_bb <- df_rho_bb$rho_true[1:nEdge]
  #truncate values outside of limits
  rho_true_bb[rho_true_bb < lims[1]] <- lims[1]
  rho_true_bb[rho_true_bb > lims[2]] <- lims[2]
  plt <- ggcorrplot2(
    values = rho_true_bb,
    colFUN = seq3div,  #user-defined function from FC_vis_funs.R
    divColor="black",
    title = paste0("Ground Truth Rho Value (", baseName, ")"),
    legTitle = 'rho value',
    lim=lims,
    diagVal=0,
    legend.position = 'right'
  )
  print(plt)
  dev.off()

  #read in reliability of ICC for T = 30 min
  load(file = file.path(dir_github, 'results', '8_BWAS', paste0('BWAS_',baseName,'_ICC_FC.RData'))) #ICC_X, ICC_X_60min
  ICC_X[ICC_X < 0] <- 0
  ICC_X_60min[ICC_X_60min < 0] <- 0
  df_rho_bb$ICC_FC_30min <- ICC_X
  df_rho_bb$ICC_FC_60min <- ICC_X_60min
  df_rho <- rbind(df_rho, df_rho_bb)

  ### Visualize ICC of FC (T = 30 min)

  ICC_X_df <- data.frame(ICC_X = ICC_X)
  ICC_X_df$ICC_group <- ifelse(ICC_X_df$ICC_X > 0.9, "excellent",
                                  ifelse(ICC_X_df$ICC_X > 0.75, "good",
                                         ifelse(ICC_X_df$ICC_X > 0.5, "moderate", "poor")))

  #calculate Q1 and Q3 of ICC_X
  quantile(ICC_X, probs = c(0.25, 0.75), na.rm = TRUE)

  pdf(file.path(dir_plots, "BWAS", paste0("ICC_FC_", baseName,".pdf")), height=4, width=4)

  Title <- paste0("ICC of FC,\nT = 30 min (", baseName, ")")

  print(ggplot(ICC_X_df, aes(x = ICC_X)) +
          geom_histogram(aes(fill = ..x..), binwidth = 0.05, color = 'gray20', breaks = seq(0,1,0.05)) +
          scale_fill_viridis_c(option = "F", direction = -1) +
          geom_vline(aes(xintercept = mean(ICC_X, na.rm = TRUE)), color = 'black', linewidth = 0.8, linetype = 'dashed') +
          scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
          cowplot::theme_cowplot() +
          xlab('ICC of FC\n(T = 30, Stringent)') +
          annotate("text", x = 0.1, y = 11000, label = paste0("Mean = ",round(mean(ICC_X, na.rm = TRUE), 2)), hjust = 0) +
          ggtitle(Title) +
          theme(legend.position = 'none',
                axis.line.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank()))

  #truncate values outside of limits
  ICC_X_lims <- ICC_X
  ICC_X_lims[ICC_X_lims < lims_ICC[1]] <- lims_ICC[1]
  ICC_X_lims[ICC_X_lims > lims_ICC[2]] <- lims_ICC[2]
  plt <- ggcorrplot2(
    values = ICC_X_lims,
    colFUN = seq7,  #user-defined function from FC_vis_funs.R
    divColor="black",
    title = Title,
    legTitle = 'ICC',
    lim=lims_ICC,
    diagVal=0,
    legend.position = 'right'
  ) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 18, margin = margin(t = 5)),
      axis.title.y = element_text(size = 18, margin = margin(t = 10)),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(x = "Node", y = "Node")

  print(plt)

  dev.off()

  ### Visualize Theoretical BWAS Proportional Bias ( = E[rho-hat]/rho)

  BWAS_propBias_math <- sqrt(ICC_y)*sqrt(ICC_X)
  print(summary(BWAS_propBias_math))
  print(quantile(BWAS_propBias_math, probs = c(0.20, 0.80), na.rm=TRUE))

  BWAS_propBias_math_mean[[baseName]] <- mean(BWAS_propBias_math, na.rm = TRUE)

  BWAS_propBias_math_df <- data.frame(BWAS_propBias_math = BWAS_propBias_math)

  pdf(file.path(dir_plots, "BWAS", paste0("BWAS_propBias_math_", baseName,".pdf")), height=4, width=4)

  Title <- paste0("Theoretical BWAS\n T = 30 min (", baseName, ")")

  print(ggplot(BWAS_propBias_math_df, aes(x = BWAS_propBias_math, fill = after_stat(x))) +
          geom_histogram(binwidth = 0.05, color = 'gray20', breaks = seq(0,1,0.05)) +
          scale_fill_viridis(option = "plasma", direction = -1, limits = c(0.5,1), name = 'rho', oob = scales::squish) +
          geom_vline(aes(xintercept = mean(BWAS_propBias_math, na.rm = TRUE)), color = 'black', linewidth = 0.8, linetype = 'dashed') +
          scale_x_continuous(labels = percent_format(accuracy = 1), breaks = seq(0,1,0.2), limits = c(0,1)) +
          scale_y_continuous(labels = percent_format(accuracy = 1)) +
          cowplot::theme_cowplot() +
          xlab('BWAS Proportional Strength\n(Total Cognition)') +
          annotate("text", x = 0.53, y = 10000, label = paste0("Mean = ",round(mean(BWAS_propBias_math, na.rm = TRUE), 2)), hjust = 0) +
          ggtitle(Title) +
          theme(legend.position = 'none',
                axis.line.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank()) +
          coord_flip())

  #truncate values outside of limits
  limits_propBias <- c(0.5, 1)
  BWAS_propBias_math_lims <- BWAS_propBias_math
  BWAS_propBias_math_lims[BWAS_propBias_math_lims < limits_propBias[1]] <- limits_propBias[1]
  BWAS_propBias_math_lims[BWAS_propBias_math_lims > limits_propBias[2]] <- limits_propBias[2]
  plt <- ggcorrplot2(
    values = BWAS_propBias_math_lims,
    colFUN = seq5_reverse,  #user-defined function from FC_vis_funs.R
    divColor="black",
    title = Title,
    legTitle = 'prop',
    lim=limits_propBias,
    diagVal=0,
    legend.position = 'right'
  ) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 18, margin = margin(t = 5)),
      axis.title.y = element_text(size = 18, margin = margin(t = 10)),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(x = "Node", y = "Node")

  print(plt)

  dev.off()

  #plot showing relationship between ICC and BWAS proportional strength

  dens <- density(ICC_X_df$ICC_X, from = 0, to = 1)
  dens_df <- data.frame(ICC_X = dens$x, dens = dens$y)
  dens_df <- dens_df %>% mutate(dens_norm = (dens - min(dens)) / (max(dens) - min(dens)))

  ICC_BWAS <- data.frame(ICC_X = ICC_X,
                         CogTotalComp = sqrt(ICC_y)*sqrt(ICC_X),
                         CogFluidComp = sqrt(ICC_demo[13])*sqrt(ICC_X),
                         ProcSpeed = sqrt(ICC_demo[23])*sqrt(ICC_X),
                         Dexterity = sqrt(ICC_demo[35])*sqrt(ICC_X))

  ICC_BWAS <- ICC_BWAS %>%
    pivot_longer(cols = -ICC_X, values_to = "BWAS_propstrength", names_to = "behavioral_measure")

  ICC_BWAS$behavioral_measure <- factor(ICC_BWAS$behavioral_measure, levels = c("CogTotalComp", "CogFluidComp", "ProcSpeed", "Dexterity"),
                                        labels = c("Total Cognition (ICC = 0.93)", "Fluid Cognition (ICC = 0.76)", "Processing Speed (ICC = 0.66)", "Dexterity (ICC = 0.49)"))

  pdf(file.path(dir_plots, "BWAS", paste0("ICC_vs_propstrength_", baseName,".pdf")), height = 4, width = 4)
  print(ggplot() +
          geom_tile(data = dens_df, aes(x = ICC_X, y = 0.5, fill = dens_norm), width = diff(dens_df$ICC_X)[1], height = 1, alpha = 0.5) +
          scale_fill_viridis_c(option = "F", direction = -1, guide = "none") +
          geom_line(data = ICC_BWAS, aes(x = ICC_X, y = BWAS_propstrength, color = behavioral_measure, group = behavioral_measure), linewidth = 0.8) +
          scale_color_viridis_d(option = "F") +
          xlab("ICC of FC") +
          ylab("BWAS Proportional Strength") +
          scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0,1,0.2), limits = c(0,1)) +
          scale_x_continuous(breaks = seq(0,1.0,0.2)) +
          cowplot::theme_cowplot() +
          theme(legend.position = c(1, 0.02),
                legend.justification = c(1, 0),
                legend.background = element_rect(fill = 'white', color = 'gray30'),
                legend.title = element_blank(),
                legend.text = element_text(size = 9.5),
                legend.margin = margin(1,4,1,4)))
  dev.off()

}

df_rho$scrub <- factor(df_rho$scrub, levels = FD_levels, labels = FD_levels_long) #None / FD5 / FD2 / FD2+
df_rho$duration <- df_rho$duration*2
df_rho$Duration2 <- paste0(df_rho$duration, " min")
df_rho$Duration2[df_rho$Duration2=="28.44 min"] <- "30 min"
df_rho$Duration2 <- factor(df_rho$Duration2, levels = paste0(seq(5,30,2.5)," min"))


################################################################

### Compute and Visualize BWAS Attenuation

df_rho$proportion <- df_rho$rho_est/df_rho$rho_true
df_rho$proportion[abs(df_rho$rho_true) < 0.01] <- NA # only consider edges with BWAS mag >= 0.1

#how many edges are included when we filter out tiny magnitude BWAS?
df_rho %>% group_by(base) %>% filter(scrub=='Stringent', duration==28.44) %>%
  summarize(mean(abs(rho_true) > 0.01))

df_rho_summ <- df_rho %>% group_by(base, scrub, duration) %>%
  summarize(proportion_mean = mean(proportion, na.rm=TRUE)) #mean over edges and visits
df_rho_summ <- filter(df_rho_summ, base!= 'P32') #exclude P32

pdf(file.path(dir_plots, "BWAS", "attenuation.pdf"), height=4, width=4.25)
for(bb in 2:nB){
  base_b <- baseNames[bb]
  p <- ggplot(filter(df_rho_summ, base==base_b),
              aes(x = duration, y = proportion_mean, group=scrub)) +
    geom_hline(yintercept = BWAS_propBias_math_mean[[base_b]], color = '#b3ce91', alpha = 0.8, linewidth = 0.8, linetype = 'dashed') +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits=c(0.25,1), breaks=seq(0.25,1,0.25)) +
    geom_line(aes(color=scrub)) +
    xlab('Scan Duration (min)') + ylab('BWAS Proportional Strength\n(Cognition Total Composite Score)') +
    scale_x_continuous(breaks=seq(10,30,10), limits = c(5,30)) +
    scale_color_manual(values = myCols, name = "Scrubbing") + scale_fill_manual(values = myCols, name = "Scrubbing") +
    cowplot::theme_cowplot() +
    theme(legend.position = c(0.95, 0.10),
          legend.justification = c("right", "bottom")) +
    ggtitle(paste0('BWAS Attenuation, (',base_b,')' ))
  print(p)
}
dev.off()


### Another version: Use bias-correct ground truth BWAS

#what is the range of correction factors (setting min value of ICC_FC to 0.1)
(maxICC_FC_36P <- max(df_rho$ICC_FC_60min[df_rho$base=='P36']))
(maxICC_FC_FIX <- max(df_rho$ICC_FC_60min[df_rho$base=='FIX']))
1 / ( sqrt(maxICC_FC_FIX) * sqrt(ICC_y) )
1 / ( sqrt(maxICC_FC_36P) * sqrt(ICC_y) )
1 / ( sqrt(0.1) * sqrt(ICC_y) )

# "ground truth" BWAS correlations with bias correction
df_rho$rho_true_corr <- df_rho$rho_true / ( sqrt(df_rho$ICC_FC_60min) * sqrt(ICC_y) )

#% of edges removed when we ignore those with ICC < 0.1
percent_removed <- df_rho %>%
  filter(duration == 28.44) %>%
  group_by(base, scrub) %>%
  summarise(
    total_edges = n(),
    removed_edges = sum(ICC_FC_60min < 0.1, na.rm = TRUE),
    percent_removed = 100 * (removed_edges / total_edges)
  ) %>%
  ungroup()

df_rho$rho_true_corr[df_rho$ICC_FC_60min < 0.1] <- NA #exclude edges with very low ICC

# get 75th, 90th, 95th, and 99th quantiles
rho_P36 <- c(unname(quantile(df_rho$rho_true_corr[df_rho$base == "P36"], probs = c(0.75, 0.90, 0.95), na.rm = TRUE)))
rho_FIX <- c(unname(quantile(df_rho$rho_true_corr[df_rho$base == "FIX"], probs = c(0.75, 0.90, 0.95), na.rm = TRUE)))
rho_FIX_GSR <- c(unname(quantile(df_rho$rho_true_corr[df_rho$base == "FIX_GSR"], probs = c(0.75, 0.90, 0.95), na.rm = TRUE)))

E_rhohat <- list()

for (bb in 2:nB) {
  baseName <- baseNames[bb]
  E_rhohat[[baseName]] <- list()
  rho_values <- if (baseName == "P36") { rho_P36 } else if (baseName == "FIX") { rho_FIX } else { rho_FIX_GSR }

  for (ss in c(1,2,5,6)) {
    scrub_ss <- FD_levels2[ss]
    load(file = file.path(dir_github, 'results', '8_BWAS', paste0('BWAS_',baseName,'_',scrub_ss,'_ICC_FC.RData')))
    ICC_X[ICC_X < 0] <- 0
    rho_edges <- outer(sqrt(ICC_X * ICC_y), rho_values, FUN = "*")
    colnames(rho_edges) <- round(rho_values, digits = 4)
    E_rhohat[[baseName]][[scrub_ss]] <- rho_edges
  }
}

#fisher z-transform of E_rhohat
z_E_rhohat <- lapply(E_rhohat, function(scrub_list) { lapply(scrub_list, psych::fisherz) })

N <- seq(500, 10000, 10)
sqrtN <- sqrt(N - 3)
rho_labels <- c("75th Percentile", "90th Percentile", "95th Percentile")

power_df <- data.frame()
for (i in seq_along(rho_labels)) {
  label <- rho_labels[i]

  for (ss in c(1,2,5,6)) {
    scrub_ss <- FD_levels2[ss]

    #get z-matrix for scrub level
    z_mat_P36 <- z_E_rhohat[["P36"]][[scrub_ss]]
    z_mat_FIX <- z_E_rhohat[["FIX"]][[scrub_ss]]
    z_mat_FIX_GSR <- z_E_rhohat[["FIX_GSR"]][[scrub_ss]]

    #expected z for each edge for each true rho value
    z_exp_P36 <- z_mat_P36[,i]
    z_exp_FIX <- z_mat_FIX[,i]
    z_exp_FIX_GSR <- z_mat_FIX_GSR[,i]

    #multiply by sqrt(N - 3) for each N
    mu_P36 <- outer(z_exp_P36, sqrtN, '*')
    mu_FIX <- outer(z_exp_FIX, sqrtN, '*')
    mu_FIX_GSR <- outer(z_exp_FIX_GSR, sqrtN, '*')

    #power = P(|Z^| * sqrt(N-3) > 1.96)
    power_P36 <- pnorm(-1.96, mu_P36, 1) + (1 - pnorm(1.96, mu_P36, 1))
    power_FIX <- pnorm(-1.96, mu_FIX, 1) + (1 - pnorm(1.96, mu_FIX, 1))
    power_FIX_GSR <- pnorm(-1.96, mu_FIX_GSR, 1) + (1 - pnorm(1.96, mu_FIX_GSR, 1))

    #summarize over edges for each N
    mean_P36 <- apply(power_P36, 2, mean, na.rm = TRUE)
    median_P36 <- apply(power_P36, 2, median, na.rm = TRUE)

    mean_FIX <- apply(power_FIX, 2, mean, na.rm = TRUE)
    median_FIX <- apply(power_FIX, 2, median, na.rm = TRUE)

    mean_FIX_GSR <- apply(power_FIX_GSR, 2, mean, na.rm = TRUE)
    median_FIX_GSR <- apply(power_FIX_GSR, 2, median, na.rm = TRUE)

    #annotate rho labels with rho value
    rho_val_P36 <- rho_P36[i]
    rho_val_FIX <- rho_FIX[i]
    rho_val_FIX_GSR <- rho_FIX_GSR[i]
    rho_label_P36 <- paste0("atop(bold('",label,"'), bold((rho == ",formatC(rho_val_P36,format="f",digits=4),")))")
    rho_label_FIX <- paste0("atop(bold('",label,"'), bold((rho == ",formatC(rho_val_FIX,format="f",digits=4),")))")
    rho_label_FIX_GSR <- paste0("atop(bold('",label,"'), bold((rho == ",formatC(rho_val_FIX_GSR,format="f",digits=4),")))")

    power_df <- rbind(power_df,
                      data.frame(N = N, mean = mean_P36, median = median_P36, 
                                 base = "P36", rho_label = rho_label_P36, scrub = scrub_ss),
                      data.frame(N = N, mean = mean_FIX, median = median_FIX, 
                                 base = "FIX", rho_label = rho_label_FIX, scrub = scrub_ss),
                      data.frame(N = N, mean = mean_FIX_GSR, median = median_FIX_GSR, 
                                 base = "FIX_GSR", rho_label = rho_label_FIX_GSR, scrub = scrub_ss))
  }
}

myCols <- c("black", "#5585bd", "#b3ce91", "#b63545") #for None / FD5 / FD2 / FD2+
power_df$scrub <- factor(power_df$scrub, levels = c("None","FD5","FD2","FD2+"), labels = FD_levels_long)

#plot sample size and power
for (bb in 2:nB) {
  baseName <- baseNames[bb]

  pdf(file.path(dir_plots, "BWAS", paste0("N_power_",baseName,".pdf")), height = 4, width = 10)
  print(ggplot(data = subset(power_df, base == baseName), aes(x = N, color = scrub, fill = scrub)) +
          geom_line(aes(y = median)) +
          facet_wrap(~ rho_label, ncol = 3, labeller = label_parsed) +
          labs(x = "Sample Size (N)", y = "Mean Power", title = "Sample Size and Power") +
          scale_color_manual(values = myCols) +
          scale_fill_manual(values = myCols) +
          scale_x_continuous(breaks = c(500,1000,2000,3000,4000,5000), limits = c(500,5000)) +
          cowplot::theme_cowplot() +
          theme(legend.position = c(0.99,0.05),
                legend.justification = c(1,0),
                legend.title = element_blank(),
                strip.background = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

for(bb in 2:nB){
  base_b <- baseNames[bb]
  png(file.path(dir_plots, "BWAS", paste0("corrected_groundTruth_",base_b,".png")), width=6, height = 5, units = "in", res = 600)
  df_rho_bb <- filter(df_rho, base==base_b)
  p <- ggplot(df_rho_bb[1:nEdge,], aes(x = rho_true, y = rho_true_corr)) +
    geom_point(aes(color = ICC_FC_60min), alpha=0.5) +
    scale_color_viridis_c('FC ICC\n(T = 60)', direction = -1) + theme_few() +
    xlab('Nominal Ground Truth rho') + ylab('Bias-Corrected Ground Truth rho') +
    ggtitle(paste0('Bias-Corrected Ground Truth BWAS (', base_b, ')'))
  print(p)
  dev.off()
}

df_rho$proportion <- df_rho$rho_est/df_rho$rho_true_corr
df_rho$proportion[abs(df_rho$rho_true_corr) < 0.01] <- NA # only consider edges with BWAS mag >= 0.1

#how many edges are included when we filter out tiny magnitude BWAS?
df_rho %>% group_by(base) %>% filter(scrub=='Stringent', duration==28.44) %>%
  summarize(mean(abs(rho_true_corr) > 0.01, na.rm=TRUE))

df_rho_summ <- df_rho %>% group_by(base, scrub, duration) %>%
  summarize(proportion_mean = mean(proportion, na.rm=TRUE)) #mean over edges and visits
df_rho_summ <- filter(df_rho_summ, base!= 'P32') #exclude P32

pdf(file.path(dir_plots, "BWAS", "attenuation_biascorr.pdf"), height=6.45, width=3)
for(bb in 2:nB){
  base_b <- baseNames[bb]
  p <- ggplot(filter(df_rho_summ, base==base_b),
              aes(x = duration, y = proportion_mean, group=scrub)) +
    geom_hline(yintercept = BWAS_propBias_math_mean[[base_b]], color = '#b3ce91', alpha = 0.8, linewidth = 0.8, linetype = 'dashed') +
    scale_y_continuous(labels = percent_format(accuracy = 1), 
                       limits=c(0.20,0.75), 
                       breaks=seq(0,75,0.1),
                       minor_breaks = seq(0,75,0.05)) +
    geom_line(aes(color=scrub)) +
    xlab('Scan Duration (min)') + ylab('BWAS Proportional Strength (mean over edges)') +
    scale_x_continuous(breaks=seq(10,30,10), limits = c(5,30)) +
    scale_color_manual(values = myCols) + scale_fill_manual(values = myCols) +
    cowplot::theme_cowplot() +
    theme(legend.position = c(0.86, 0.10),
          legend.justification = c("right", "bottom"),
          legend.title = element_blank(),
          panel.grid.minor.y = element_line(color = "gray80", size = 0.20),
          panel.grid.major.y = element_line(color = "gray80", size = 0.20)) +
    ggtitle(paste0('BWAS Attenuation, Bias-Corrected GT (',base_b,')' ))

  if (base_b == "P36") {
    p <- p + annotate("text", x = 6, y = BWAS_propBias_math_mean[[base_b]] + 0.030, #+ 0.015,
                      label = "Theoretical (T = 30)", color = '#b3ce91', hjust = 0)
  }

  print(p)
}
dev.off()
