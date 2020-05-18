# TODO(AW) make this reproducible!

## a script to investigate the coverage of my confidence set procedure
## HIV dataset

##### first: run bash script

### HIV dataset
setwd("Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/")
source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/miscellaneous_phylogenetics_programs.R")

################################################
#################### n=20
################################################
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated20/split/trees/")
length(trees) ## should be 20000
for (tree in trees) {
  path <- paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated20/split/trees/", tree, sep="")
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated20/trees/", tree,sep=""), row.names = F, quote = F, col.names = F)
}
### create 1000 batches of 20 trees, run the program on each batch, see if covers
batches <- 1000
containsTruth <- rep(NA, batches)
containsTruthAbs <- rep(NA, batches); mu0Abs <- c(0.002854010, 0.0379729)
batchSize <- length(trees)/ batches
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated20/trees/", full.names = T)
pts <- matrix(NA, nrow=batches, ncol=2)
pts_true <- matrix(NA, nrow=batches, ncol=2)
covs <- list()
for (batchNumber in 1:batches) {
  ### create the batch
  path <- trees[(batchNumber - 1)*batchSize + 1]
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/test1.txt"), row.names = F, quote = F, col.names = F, append = F)
  for (j in 2:batchSize) {
    path <- trees[(batchNumber - 1)*batchSize + j]
    tree1 <- as.character(read.table(path)[1,1])
    write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/test1.txt"), row.names = F, quote = F, col.names = F, append = T)
  }

  ### run the script to get the log maps, and load them
  system('java -jar coverage_java.jar')
  source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/maps.R")

  ### get the log map of the true tree
  system('java -jar coverage_truth_java.jar')
  source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/map_truth.R")

  ### check if true tree is in confidence set
  observations <- matrix(NA, ncol = length(tree1coveragesim), nrow = batchSize)
  for (i in 1:batchSize) {
    observations[i, ] <- get(paste("tree", i, "coveragesim", sep=""))
  }
  mu0 <- trueTreeLogMap
  xbar = apply(observations, 2,  mean)
  S = cov(observations)
  n=batchSize; m=length(tree1coveragesim)
  containsTruth[batchNumber] <- t(xbar - mu0)%*%solve(S)%*%(xbar - mu0)
  containsTruthAbs[batchNumber] <- t(xbar - mu0Abs)%*%solve(S)%*%(xbar - mu0Abs)
  # pts[batchNumber, ] <- xbar
  # pts_true[batchNumber, ] <- mu0
  # covs[[batchNumber]] <- S
}
n=batchSize; m=length(tree1coveragesim)
c(mean(containsTruth < (m*(n-1))/(n*(n-m))*qf(0.90, m, n-m)),
  mean(containsTruth < (m*(n-1))/(n*(n-m))*qf(0.95, m, n-m)),
  mean(containsTruth < (m*(n-1))/(n*(n-m))*qf(0.99, m, n-m)),
  mean(containsTruth < (m*(n-1))/(n*(n-m))*qf(0.999, m, n-m)))
containsTruth_HIV_20 <- containsTruth

################################################
#################### n=50
################################################
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated50/split/trees/")
length(trees) ## should be 50000
for (tree in trees) {
  path <- paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated50/split/trees/", tree, sep="")
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated50/trees/", tree,sep=""), row.names = F, quote = F, col.names = F)
}
batches <- 1000
containsTruth <- rep(NA, batches)
containsTruthAbs <- rep(NA, batches); mu0Abs <- c(0.002854010, 0.0379729)
batchSize <- length(trees)/ batches
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated50/trees/", full.names = T)
for (batchNumber in 1:batches) {
  path <- trees[(batchNumber - 1)*batchSize + 1]
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/test1.txt"), row.names = F, quote = F, col.names = F, append = F)
  for (j in 2:batchSize) {
    path <- trees[(batchNumber - 1)*batchSize + j]
    tree1 <- as.character(read.table(path)[1,1])
    write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/test1.txt"), row.names = F, quote = F, col.names = F, append = T)
  }
  system('java -jar coverage_java.jar')
  source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/maps.R")
  system('java -jar coverage_truth_java.jar')
  source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/map_truth.R")
  observations <- matrix(NA, ncol = length(tree1coveragesim), nrow = batchSize)
  for (i in 1:batchSize) {
    observations[i, ] <- get(paste("tree", i, "coveragesim", sep=""))
  }
  mu0 <- trueTreeLogMap
  xbar = apply(observations, 2,  mean)
  S = cov(observations)
  n=batchSize; m=length(tree1coveragesim)
  containsTruth[batchNumber] <- t(xbar - mu0)%*%solve(S)%*%(xbar - mu0)
  containsTruthAbs[batchNumber] <- t(xbar - mu0Abs)%*%solve(S)%*%(xbar - mu0Abs)
}
n=batchSize; m=length(tree1coveragesim)
c(mean(containsTruth < (m*(n-1))/(n*(n-m))*qf(0.90, m, n-m)),
  mean(containsTruth < (m*(n-1))/(n*(n-m))*qf(0.95, m, n-m)),
  mean(containsTruth < (m*(n-1))/(n*(n-m))*qf(0.99, m, n-m)),
  mean(containsTruth < (m*(n-1))/(n*(n-m))*qf(0.999, m, n-m)))
containsTruth_HIV_50 <- containsTruth

################################################
#################### n=100
################################################
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated100/split/trees/")
length(trees) ## should be 100000
for (tree in trees) {
  path <- paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated100/split/trees/", tree, sep="")
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated100/trees/", tree,sep=""), row.names = F, quote = F, col.names = F)
}
batches <- 1000
containsTruth <- rep(NA, batches)
containsTruthAbs <- rep(NA, batches); mu0Abs <- c(0.002854010, 0.0379729)
batchSize <- length(trees)/ batches
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_hiv/generated100/trees/", full.names = T)
for (batchNumber in 1:batches) {
  ## skip 725 & 815, not working?
  path <- trees[(batchNumber - 1)*batchSize + 1]
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/test1.txt"), row.names = F, quote = F, col.names = F, append = F)
  for (j in 2:batchSize) {
    path <- trees[(batchNumber - 1)*batchSize + j]
    tree1 <- as.character(read.table(path)[1,1])
    write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/test1.txt"), row.names = F, quote = F, col.names = F, append = T)
  }
  system('java -jar coverage_java.jar')
  source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/maps.R")
  system('java -jar coverage_truth_java.jar')
  source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/map_truth.R")
  observations <- matrix(NA, ncol = length(tree1coveragesim), nrow = batchSize)
  for (i in 1:batchSize) {
    observations[i, ] <- get(paste("tree", i, "coveragesim", sep=""))
  }
  mu0 <- trueTreeLogMap
  xbar = apply(observations, 2,  mean)
  S = cov(observations)
  n=batchSize; m=length(tree1coveragesim)
  containsTruth[batchNumber] <- t(xbar - mu0)%*%solve(S)%*%(xbar - mu0)
  containsTruthAbs[batchNumber] <- t(xbar - mu0Abs)%*%solve(S)%*%(xbar - mu0Abs)
}
n=batchSize; m=length(tree1coveragesim)

containsTruth_HIV_100 <- containsTruth
m=2
n=20
100*c(mean(containsTruth_HIV_20 < (m*(n-1))/(n*(n-m))*qf(0.90, m, n-m)),
  mean(containsTruth_HIV_20 < (m*(n-1))/(n*(n-m))*qf(0.95, m, n-m)),
  mean(containsTruth_HIV_20 < (m*(n-1))/(n*(n-m))*qf(0.99, m, n-m)),
  mean(containsTruth_HIV_20 < (m*(n-1))/(n*(n-m))*qf(0.999, m, n-m)))
n=50
100*c(mean(containsTruth_HIV_50 < (m*(n-1))/(n*(n-m))*qf(0.90, m, n-m)),
  mean(containsTruth_HIV_50 < (m*(n-1))/(n*(n-m))*qf(0.95, m, n-m)),
  mean(containsTruth_HIV_50 < (m*(n-1))/(n*(n-m))*qf(0.99, m, n-m)),
  mean(containsTruth_HIV_50 < (m*(n-1))/(n*(n-m))*qf(0.999, m, n-m)))
n=100
100*c(mean(containsTruth_HIV_100 < (m*(n-1))/(n*(n-m))*qf(0.90, m, n-m)),
  mean(containsTruth_HIV_100 < (m*(n-1))/(n*(n-m))*qf(0.95, m, n-m)),
  mean(containsTruth_HIV_100 < (m*(n-1))/(n*(n-m))*qf(0.99, m, n-m)),
  mean(containsTruth_HIV_100 < (m*(n-1))/(n*(n-m))*qf(0.999, m, n-m)))


### Old stuff from debugging
## Soln turned out to be brackets error. Doh!


# ### sim 1
# lim1 <- 0.5
# plot(0,0,xlim = c(-lim1, lim1), ylim=c(-lim1, lim1))
# points(observations)
# points(xbar[1], xbar[2], pch=16)
# points(mu0[1], mu0[2], col="red")
# grid1 <- cbind(rep(seq(min(observations[,1]), max(observations[,1]), length=100), each=100),
#                rep(seq(min(observations[,2]), max(observations[,2]), length=100), 100))
# grid1 <- cbind(rep(seq(-lim1, lim1, length=100), each=100),
#                rep(seq(-lim1, lim1, length=100), 100))
# for (i in 1:dim(grid1)[1]) {
#   tmp <- t(xbar - grid1[i, ])%*%solve(S)%*%(xbar - grid1[i, ])
#   if (tmp < (m*(n-1))*qf(0.95, m, n-m)/n*(n-m)) points(grid1[i, 1], grid1[i, 2], col="grey")
#   if (tmp < qf(0.95, m, n-m)) points(grid1[i, 1], grid1[i, 2], col="blue")
# }
# ## conf int way too large
# sqrt(m*(n-1)*qf(alpha, m, n-m)/(n*(n-m)))
#
# dim((xbar - grid1)%*%solve(S)%*%t(xbar - grid1))
#
# points(grid1)
# ###
# plot(pts, cex=0.5, pch=16)
# points(pts_true, cex=0.5, pch=16, col="blue")
# points(0.002854010, 0.0379729, col="red",pch=16) ## true mean
# points(apply(pts, 2, mean)[1],apply(pts, 2, mean)[2], col="green", pch=16)
