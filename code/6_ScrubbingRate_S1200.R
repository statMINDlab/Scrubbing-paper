####################################################################
# Who: Damon Pham
# When: Spring 2024
# What: Determine the number of subjects with at least 30 min after FD2 scrubbing
# How: First, 3_AggFD.R must be run
# Where: Run on HPC, since it requires access to cluster storage files
####################################################################

## ---------------------------------------------------------------------------
source("0_SharedCode.R")
setwd(dir_slate)

x <- readRDS(file.path(dir_results, "results/3_AggFD/withS1200_AggFD.rds"))

## Compute the number of time points exceeding modFD = 0.2mm (Note: this does NOT include DVARS)
x1 <- x
x1[is.na(x1)] <- Inf
z <- (1185 - rowSums(x1 > .2)) * .72 / 60
z <- colSums(matrix(z, nrow=4)) #average over sessions for a given subject
hist(z, breaks=30)
table(z > 30) #how many subjects have at least 30 min remaining after FD2 scrubbing?
