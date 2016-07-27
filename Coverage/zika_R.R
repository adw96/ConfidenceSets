## a script to investigate the coverage of my confidence set procedure
## zika dataset

##### first: run bash script

### zika dataset
setwd("Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/")
source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/miscellaneous_phylogenetics_programs.R")

################################################
#################### n=20 
################################################
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated20/split/trees/")
length(trees) ## should be 20000
for (tree in trees) {
  path <- paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated20/split/trees/", tree, sep="")
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated20/trees/", tree,sep=""), row.names = F, quote = F, col.names = F)
}
### create 1000 batches of 20 trees, run the program on each batch, see if covers
batches <- 1000 
containsTruth <- rep(NA, batches)
batchSize <- length(trees)/ batches
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated20/trees", full.names = T)
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
  system('java -jar coverage_truth_zika_java.jar')
  source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/map_truth.R")
  
  ### check if true tree is in confidence set
  observations <- matrix(NA, ncol = length(tree1coveragesim), nrow = batchSize)
  for (i in 1:batchSize) {
    observations[i, ] <- get(paste("tree", i, "coveragesim", sep=""))
  }
  mu0 <- trueTreeLogMap
  xbar = apply(observations, 2,  mean)
  S = cov(observations)
  containsTruth[batchNumber] <- t(xbar - mu0)%*%solve(S)%*%(xbar - mu0)
}
containsTruth_Zika_20 <- containsTruth


################################################
#################### n=50 
################################################
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated50/split/trees/")
length(trees) ## should be 50000
for (tree in trees) {
  path <- paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated50/split/trees/", tree, sep="")
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated50/trees/", tree,sep=""), row.names = F, quote = F, col.names = F)
}
batches <- 1000
containsTruth <- rep(NA, batches)
batchSize <- length(trees)/ batches
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated50/trees/", full.names = T)
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
  system('java -jar coverage_truth_zika_java.jar')
  source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/map_truth.R")
  observations <- matrix(NA, ncol = length(tree1coveragesim), nrow = batchSize)
  for (i in 1:batchSize) {
    observations[i, ] <- get(paste("tree", i, "coveragesim", sep=""))
  }
  mu0 <- trueTreeLogMap
  xbar = apply(observations, 2,  mean)
  S = cov(observations)
  containsTruth[batchNumber] <- t(xbar - mu0)%*%solve(S)%*%(xbar - mu0) 
}
containsTruth_Zika_50 <- containsTruth

################################################
#################### n=100 
################################################
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated100/split/trees/")
length(trees) ## should be 100000
for (tree in trees) {
  path <- paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated100/split/trees/", tree, sep="")
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated100/trees/", tree,sep=""), row.names = F, quote = F, col.names = F)
}
batches <- 1000
containsTruth <- rep(NA, batches)
batchSize <- length(trees)/ batches
trees <- list.files(path = "/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/coverage_zika/generated100/trees/", full.names = T)
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
  system('java -jar coverage_truth_zika_java.jar')
  source("/Users/adw96/Documents/Phylogenetics/Tree Inference/Tree Inference Analysis/Coverage/generated/map_truth.R")
  observations <- matrix(NA, ncol = length(tree1coveragesim), nrow = batchSize)
  for (i in 1:batchSize) {
    observations[i, ] <- get(paste("tree", i, "coveragesim", sep=""))
  }
  mu0 <- trueTreeLogMap
  xbar = apply(observations, 2,  mean)
  S = cov(observations)
  containsTruth[batchNumber] <- t(xbar - mu0)%*%solve(S)%*%(xbar - mu0) 
}
containsTruth_Zika_100 <- containsTruth

m=length(tree1coveragesim)
n=20
100*c(mean(containsTruth_Zika_20 < (m*(n-1))/(n*(n-m))*qf(0.90, m, n-m)),
      mean(containsTruth_Zika_20 < (m*(n-1))/(n*(n-m))*qf(0.95, m, n-m)), 
      mean(containsTruth_Zika_20 < (m*(n-1))/(n*(n-m))*qf(0.99, m, n-m)),
      mean(containsTruth_Zika_20 < (m*(n-1))/(n*(n-m))*qf(0.999, m, n-m)))
n=50
100*c(mean(containsTruth_Zika_50 < (m*(n-1))/(n*(n-m))*qf(0.90, m, n-m)),
      mean(containsTruth_Zika_50 < (m*(n-1))/(n*(n-m))*qf(0.95, m, n-m)), 
      mean(containsTruth_Zika_50 < (m*(n-1))/(n*(n-m))*qf(0.99, m, n-m)),
      mean(containsTruth_Zika_50 < (m*(n-1))/(n*(n-m))*qf(0.999, m, n-m)))
n=100
100*c(mean(containsTruth_Zika_100 < (m*(n-1))/(n*(n-m))*qf(0.90, m, n-m)),
      mean(containsTruth_Zika_100 < (m*(n-1))/(n*(n-m))*qf(0.95, m, n-m)), 
      mean(containsTruth_Zika_100 < (m*(n-1))/(n*(n-m))*qf(0.99, m, n-m)),
      mean(containsTruth_Zika_100 < (m*(n-1))/(n*(n-m))*qf(0.999, m, n-m)))
