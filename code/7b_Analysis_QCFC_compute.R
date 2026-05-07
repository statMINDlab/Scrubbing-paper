####################################################################
# Who: Amanda Mejia
# When: Fall 2024
# What: Compute QC-FC (naive and partitioned)
# How: First, 3_AggFD.R and 4_AggFCandFlags.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## ---------------------------------------------------------------------------
source("0_SharedCode.R")
setwd(dir_slate)


### COMPUTE QC-FC -------------------------------------------------------------

### 1. Get mean FD for each session

#full HCP
x <- readRDS(file.path(dir_slate, "results/3_AggFD/withS1200_AggFD.rds"))
meanFD <- rowMeans(x, na.rm = FALSE) #sessions with any NA values will be set to NA (affects 240 sessions, 217 of which are NA for all time points)
meanFD <- matrix(meanFD, nrow = 4) #runs x subjects

#remove sessions with less than 1 minute of unflagged scan duration (~83 volumes) for FD2+
flags <- readRDS(file.path(dir_results, 'results/4_AggFC/withS1200_withDVARS_flags_P36_FD___og_nfc_l4___0.2_plus.rds'))
flags_14.22 <- flags[,,,"first_14.22_mins",]
unflagged_count <- apply(!flags_14.22, c(1,2), sum, na.rm = TRUE) #runs x subjects
remove <- which(unflagged_count <= 83, arr.ind = TRUE) # num excluded sessions x 2 (row/col positions in )

meanFD[remove] <- NA #set these sessions to NA

#how many subjects missing sessions?
table(colSums(is.na(meanFD)))
# 0   1   2   3   4
# 938  52  86  11  26
# before excluding sessions with less than 1 min:
#    0    1    2    3    4
# 1004   20   65    6   18

#get subject names
tmp <- readRDS(file.path(dir_slate, 'results', '4_AggFC', "withS1200_withDVARS_FC_FIX_Base_first14.22mins.rds"))
subjects <- dimnames(tmp)[[2]]
runs <- dimnames(tmp)[[1]] #"test_LR1" "test_LR2" "test_RL1" "test_RL2"

meanFD_df <- as.data.frame(t(meanFD), row.names = subjects)
names(meanFD_df) <- runs
meanFD_df$subject <- subjects
meanFD_df <- reshape2::melt(meanFD_df, id.vars='subject', variable.name = 'run', value.name = 'meanFD')

#for naive within-subject QC-FC:

#a) compute mean and variance of meanFD ACROSS subjects
meanFD_all <- colMeans(meanFD, na.rm=TRUE) #mean of each subject, across visits
meanFD_var <- var(meanFD_all, na.rm=TRUE) #var over subjects

#b) group subjects by their average motion level
# meanFD_qrt <- quantile(meanFD_all, c(0.25,0.5,0.75), na.rm=TRUE)
# hist(meanFD_all, breaks=30)
# abline(v = meanFD_qrt, col='red', lwd=2)
# meanFD_grp <- ifelse(meanFD_all < meanFD_qrt[1], 1,
#                      ifelse(meanFD_all < meanFD_qrt[2], 2,
#                             ifelse(meanFD_all < meanFD_qrt[3], 3, 4)))
# meanFD_hiFD <- meanFD_grp %in% 3:4
# meanFD_loFD <- meanFD_grp %in% 1:2

#c) sanity check: compute mean and variance of meanFD WITHIN subjects (use most extreme sessions from each subject)
meanFD_min <- apply(meanFD, 2, min, na.rm=TRUE); meanFD_min[is.infinite(meanFD_min)] <- NA #this happens when a subject is missing
meanFD_max <- apply(meanFD, 2, max, na.rm=TRUE); meanFD_max[is.infinite(meanFD_max)] <- NA #this happens when a subject is missing
meanFD_change <- meanFD_max - meanFD_min
meanFD_ismin <- (meanFD == matrix(rep(meanFD_min, each = 4), nrow=4))
meanFD_ismax <- (meanFD == matrix(rep(meanFD_max, each = 4), nrow=4))
meanFD_var_within <- 0.5*(meanFD_max - meanFD_min)^2 #var of each subject, across most extreme visits
meanFD_var_within0 <- apply(meanFD, 2, var, na.rm=TRUE) #var of each subject, across all visits

#d) identify subjects with high variance of motion (similar to population variance)
plot(meanFD_var_within0, meanFD_var_within,
     xlab='Within-subject Variance across all Visits',
     ylab='Within-subject Variance across extreme Visits',
     xlim=c(0,0.2), ylim=c(0,0.2))
abline(a=0,b=1,lty=2)
abline(h = meanFD_var, col='red') #variance across subjects

meanFDvar0_qrt <- quantile(meanFD_var_within0, c(0.25,0.5,0.75), na.rm=TRUE)
meanFDvar_qrt <- quantile(meanFD_var_within, c(0.25,0.5,0.75), na.rm=TRUE)

#takeaway: the subjects with above-median meanFD_var_within (var across extreme sessions) are very similar to population-level variance (their median is close to the population var)

pdf(file.path(dir_github, "plots", "QCFC","var_of_meanFD.pdf"), width=6, height=8)
par(mfrow=c(2,1))
hist(sqrt(meanFD_var_within), breaks=100, xlim=c(0,0.6),
     main = "Variance across Extreme Sessions", xlab = "Within-Subject SD")
abline(v = sqrt(meanFD_var), col='red') #pop-level var
abline(v = sqrt(meanFDvar_qrt[2]), col='royalblue1', lty=2)
abline(v = sqrt(meanFDvar_qrt[3]), col='royalblue1', lty=1)
legend('topright', legend = c('Population SD','Median of Within-Subject SD','Q3'), col=c('red','royalblue1','royalblue1'), lty=c(1,2,1))
hist(sqrt(meanFD_var_within0), breaks=100, xlim=c(0,0.6),
     main = "Variance across All Sessions", xlab = "Within-Subject SD")
abline(v = sqrt(meanFD_var), col='red')
abline(v = sqrt(meanFDvar0_qrt[2]), col='royalblue1', lty=2)
abline(v = sqrt(meanFDvar0_qrt[3]), col='royalblue1', lty=1)
legend('topright', legend = c('Population SD','Median of Within-Subject SD','Q3'), col=c('red','royalblue1','royalblue1'), lty=c(1,2,1))
dev.off()

#identify high within-subject variance subjects, for a sanity check of partitioned QC-FC
meanFD_hivar <- (meanFD_var_within > quantile(meanFD_var_within, 0.65, na.rm=TRUE))

pdf(file.path(dir_github, "plots", "QCFC","var_of_meanFD2.pdf"), width=6, height=8)
par(mfrow=c(2,1))
meanFD_diff <- meanFD_max - meanFD_min
hist((meanFD_all), breaks=(seq(0,0.9,0.01)), freq = FALSE, ylim=c(0,15), col=alpha('black', 0.5), main = 'Distribution of Mean FD', xlab = 'Mean FD', xlim=c(0,0.9))
hist((meanFD_diff), breaks=(seq(0,0.9,0.01)), add=TRUE, col=alpha('hotpink', 0.5), freq = FALSE)
hist((meanFD_diff[meanFD_hivar]), breaks=(seq(0,0.9,0.01)), add=TRUE, col=alpha('blue', 0.5), freq = FALSE)
legend('topright', legend = c('population distribution','within-subject differences', 'within-subject differences (hivar subjects)'),
       col = alpha(c('black','hotpink','blue'), 0.5), lwd=5)
y <- c(meanFD_all, meanFD_diff[meanFD_hivar])
grp <- c(rep('pop', length(meanFD_all)), rep('within (hivar)', length(meanFD_diff[meanFD_hivar])))
boxplot(y ~ grp, col = alpha(c('black','blue'), 0.5), horizontal=TRUE, xlab = 'Mean FD', ylim=c(0,0.9))
dev.off()

### 2. Partition FD into between-subject and within-subject variables

FD_vec <- c(meanFD)
FD_vec_std <- scale(FD_vec) #normalize to mean 0, var 1

#for partitioned QC-FC: separate the within-subject and between-subject effects
#between-subject effects
x_i_bar <- colMeans(meanFD, na.rm=TRUE) #mean of each subject, across visits
x_bar <- mean(meanFD, na.rm=TRUE) #overall mean
x_between <- x_i_bar - x_bar
#within-subject effects
x_within <- scale(meanFD, scale=FALSE) #subtract the within-subject mean (among non-NA visits)
x_within_vec <- c(x_within) #x_within is nSess x nSubj , so x_within_vec grouped by subject
x_between_vec <- rep(x_between, each=4) #repeat each subject's value nSess times
#form design matrix, center and scale for SD=1
X <- cbind(scale(x_within_vec), scale(x_between_vec))

#compute correction factor for repeated measures QC-FC to obtain partial correlations
corXZ <- cor(x_between_vec, x_within_vec, use = "complete.obs") #appears on the numerator of the correction factor; should be close to zero

#we will Fisher z-transform the FC values prior to computing QC-FC

#QCFC_naive_df --> cor(FD, FC)
#QCFC_adj_df --> partition QC-FC correlations into within- and between-subject associations using multiple linear regression
#QCFC_within_df --> compute within-subject associations in a simpler way as a sanity check

### 3. Get FC estimates (full duration) by scrubbing/base method

#loop over base denoising methods
for (bb in 1:nB) {
  baseName <- baseNames[bb]
  cat(baseName, "\n\n")

  QCFC_naive_df <- QCFC_adj_df <- QCFC_within_df <- NULL

  #loop over scrubbing levels (None, FD5, FD2, FD2+)
  #for (ss in c(1,2,5,6)) {
  for (ss in 1:6) {
    pfname_ss <- scrubNames[ss] #partial file name
    scrub_ss <- FD_levels2[ss]
    cat(scrub_ss, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    fname_bs <- paste0('withS1200_withDVARS_FC_', baseName, '_', pfname_ss, '.rds')
    FC_bs <- readRDS(file.path(dir_results, 'results', '4_AggFC', fname_bs))[,,,1,1] #sess x sub x edge

    #excluding runs with < 1 min remaining for FD2+

    for (s in unique(remove[,1])) { #unique(remove[,1]) indexes runs (i.e. retest_LR1)
      subjects_remove_s <- remove[remove[,1] == s, 2] #subjects to remove for run type s
      FC_bs[s, subjects_remove_s, ] <- NA
    }

    ### NAIVE

    # #for each edge, compute cor(FD, FC) across all sessions the traditional way
    # QCFC0_bs <- apply(FC_bs, 3, function(y){
    #   FC_vec <- fishZ(c(y)) #vectorize 4xnSess matrix
    #   cor(FC_vec, FD_vec, use = 'complete')
    # })

    #compute cor(FD, FC) using linear regression (this is EQUAL to cor(x,y) since SD(x)=SD(y)=1)
    QCFC_bs <- apply(FC_bs, 3, function(y){
      FC_vec <- fishZ(c(y)) #vectorize 4xnSess matrix, vector is grouped by subject
      FC_vec_std <- scale(FC_vec) #normalize to mean 0, var 1
      mod <- lm(FC_vec_std ~ FD_vec_std - 1)
      c(mod$coefficients[1], #coefficient estimate = Pearson correlation
        summary(mod)$coefficients[,4]) #p-value
    })
    QCFC_naive_df_bs <- data.frame(QCFC = QCFC_bs[1,],
                                   pval = QCFC_bs[2,],
                                   base = baseName,
                                   scrub = scrub_ss)
    QCFC_naive_df <- rbind(QCFC_naive_df, QCFC_naive_df_bs)
    #all.equal(QCFC0_bs, QCFC_bs) #TRUE (pearson correlation = SLR regression coefficient)

    ### PARTITIONED

    #helper function to correct MLR coefficients to obtain partial correlations
    MLR2pcor <- function(coefX, corXZ, Y, Z){
      #Y represents FC
      #X represents FD_predictor (this is the FD variable currently of focus, either within- or between-subject)
      #Z represents FD_nuisance (this is the other FD variable)
      #coefX is the MLR coefficient for FD_predictor
      #corXZ = cor(FD_between, FD_within) (computed at the top, same for all edges)
      corYZ <- cor(Y, Z, use = "complete.obs") #cor(FC, FD_nuisance) (edge-specific)
      valnum <- 1 - corXZ^2
      valdenom <- 1 - corYZ^2
      correction <- sqrt(valnum/valdenom)
      return(coefX * correction)
    }

    #partitioned QC-FC: separate the within-subject and between-subject effects
    QCFC_adj_bs <- apply(FC_bs, 3, function(y){
      FC_vec <- fishZ(c(y)) #vectorize 4xnSubj matrix, vector is grouped by subject
      FC_vec_std <- scale(FC_vec) #normalize to mean 0, var 1

      #obtain MLR coefficients
      mod <- lm(FC_vec_std ~ X - 1)
      coef_MLR <- coefficients(mod)

      #perform correction to obtain partial correlations
      pcor_X1 <- MLR2pcor(coef_MLR[1], corXZ, Y = FC_vec_std, Z = X[,2]) #within-subject
      pcor_X2 <- MLR2pcor(coef_MLR[2], corXZ, Y = FC_vec_std, Z = X[,1]) #between-subject

      result <- c(coef_MLR, #within-subject, between-subject
        pcor_X1, pcor_X2, #within-subject, between-subject
        summary(mod)$coefficients[,4]) #p-values (same for coefficients and partial correlations)
      names(result) <- c('MLR_within', 'MLR_between', 'pcor_within', 'pcor_between', 'pval_within', 'pval_between')
      return(result)

    })

    QCFC_adj_df_bs <- data.frame(QCFC0_within = QCFC_adj_bs[1,], QCFC0_between = QCFC_adj_bs[2,], #uncorrected effects from MLR
                                 QCFC_within = QCFC_adj_bs[3,], QCFC_between = QCFC_adj_bs[4,], #partial correlations after correction
                                 pval_within = QCFC_adj_bs[5,], pval_between = QCFC_adj_bs[6,],
                                 base = baseName,
                                 scrub = scrub_ss)
    QCFC_adj_df <- rbind(QCFC_adj_df, QCFC_adj_df_bs)

    ### within-subject QC-FC: compute just within-subject effect in a simple way

    FC_bs_minFD <- apply(FC_bs, 3, function(x) colSums(x*meanFD_ismin)) #set all FC values to zero except the one coinciding with the min-motion session
    FC_bs_maxFD <- apply(FC_bs, 3, function(x) colSums(x*meanFD_ismax)) #set all FC values to zero except the one coinciding with the max-motion session
    FC_bs_change <- fishZ(FC_bs_maxFD) - fishZ(FC_bs_minFD)
    QCFC_within_bs <- apply(FC_bs_change, 2, function(x){
      cor_all <- cor(x, meanFD_change, use = 'complete')
      cor_hivar <- cor(x[meanFD_hivar], meanFD_change[meanFD_hivar], use = 'complete')
      #cor_hivar_hiFD <- cor(x[meanFD_hivar & meanFD_hiFD], meanFD_change[meanFD_hivar & meanFD_hiFD], use = 'complete')
      #cor_hivar_loFD <- cor(x[meanFD_hivar & meanFD_loFD], meanFD_change[meanFD_hivar & meanFD_loFD], use = 'complete')
      return(c(cor_all, cor_hivar)) #, cor_hivar_hiFD, cor_hivar_loFD))
    })
    QCFC_within_df_bs <- as.data.frame(t(QCFC_within_bs))
    names(QCFC_within_df_bs) <- c('all','hivar') #,'hivar_hiFD','hivar_loFD')
    QCFC_within_df_bs$base <- baseName
    QCFC_within_df_bs$scrub <- scrub_ss
    QCFC_within_df <- rbind(QCFC_within_df, QCFC_within_df_bs)

  } #end loop over scrubbing levels

  save(QCFC_naive_df, QCFC_adj_df, QCFC_within_df,
       file = file.path(dir_github, 'results', '7_QCFC', paste0('QCFC_df_',baseName,'.RData')))

} #end loop over base denoising methods

