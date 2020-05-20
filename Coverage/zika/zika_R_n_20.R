## A script to investigate the coverage of my confidence set procedure
## for the Zika dataset for n=20

## Amy Willis, 2016

## Update May 2020: This script was absolutely irreproducible. I have modified it
## to use local paths and only distributed programs. I have not rerun it
## to confirm that the output is identical, but I hope that the increase in usability
## will allow others to confirm or deny. Please do not hesitate to reach out
## with questions, requests or comments via https://github.com/adw96/ConfidenceSets/issues.

## You should be able to clone the repo github.com/adw96/ConfidenceSets/
## and run this script in an R project in the subdirectory Coverage/zika.

# To do this analysis for n=50 and n=100, run this script with all instances of
# n=20 replaced with the appropriate n.

source("../miscellaneous_phylogenetics_programs.R")

# Step 0: run bash script `zika_shell.sh` to simulate trees

# Step 1: Pull out only the the tree output of PhyML into new files

trees <- list.files(path = "generated20/split/trees/")
length(trees) ## should be 20000
for (tree in trees) {
  path <- paste("generated20/split/trees/", tree, sep="")
  tree1 <- as.character(read.table(path)[1,1])
  write.table(convert_tree(tree1), paste("generated20/trees/", tree,sep=""),
              row.names = F, quote = F, col.names = F)
}

## Step 2: run the Bacak algorithm for lots and lots of iterations (>10^6) on these trees
## and custom check for convergence. Save resulting true tree as base_tree_zika.txt.


### Step 3: find m, call it mm
system('java -jar ../../TreePrograms/logmap.jar base_tree_zika.txt base_tree_zika.txt | tail -n 1 > base_tree_log_map.R')
source("base_tree_log_map.R")
mm <- length(logMap)


## Step 4: create 1000 batches of 20 trees. Find the LM from
## base tree to new tree. Calculate F score for each set of 20 observations.
batches <- 1000
containsTruth <- rep(NA, batches)
batchSize <- length(trees)/ batches
trees <- list.files(path = "generated20/trees", full.names = T)
for (batchNumber in 1:100) {

  
  observations <- matrix(NA, ncol = mm, nrow = batchSize)
  
  # write the trees in the batch to a file 
  path <- trees[(batchNumber - 1)*batchSize + 1]
  the_tree <- as.character(read.table(path)[1,1])
  write.table(convert_tree(the_tree),
              paste("the_sample_of_trees_tmp.txt"),
              row.names = F, quote = F, col.names = F, append = F)
  for (j in 2:batchSize) {
    # write the tree to output
    path <- trees[(batchNumber - 1)*batchSize + j]
    the_tree <- as.character(read.table(path)[1,1])
    write.table(convert_tree(the_tree),
                paste("the_sample_of_trees_tmp.txt"),
                row.names = F, quote = F, col.names = F, append = T)

  }
  
  # find the sample Frechet mean
  system('java -jar ../../TreePrograms/meantree.jar the_sample_of_trees_tmp.txt 0.0001 1000 > sample_mean_tmp.R')
  system('cat sample_mean_tmp.R | head -n 1 > sample_mean_tree_tmp.txt')
  system('cat sample_mean_tmp.R | tail -n 1 > sample_mean_lr_tmp.R')
  
  source("sample_mean_lr_tmp.R")
  tree0 # this is the log map of the sample mean
  
  # find log map of the base tree from the sample mean
  system('java -jar ../../TreePrograms/logmap.jar sample_mean_tree_tmp.txt base_tree_zika.txt | tail -n 1 > new_tree_log_map.R')
  source("new_tree_log_map.R")
  mu0 <- logMap # log map of true tree
   
  # find log map of each tree in the sample from the sample mean
  for (j in 1:batchSize) {
    system(paste("cat the_sample_of_trees_tmp.txt | head -n", j, " | tail -n 1 > jth_tree_tmp.txt"))
    system('java -jar ../../TreePrograms/logmap.jar sample_mean_tree_tmp.txt jth_tree_tmp.txt | tail -n 1 > new_tree_log_map.R')
    source("new_tree_log_map.R")
    logMap
    observations[j, ] <- logMap
  }
  observations
  
  xbar = apply(observations, 2,  mean)
  S = cov(observations)
  containsTruth[batchNumber] <- t(xbar - mu0)%*%solve(S)%*%(xbar - mu0)
}

nn=20
100*c(mean(containsTruth < (mm*(nn-1))/(n*(nn-mm))*qf(0.90, mm, nn-mm)),
      mean(containsTruth < (mm*(nn-1))/(n*(nn-mm))*qf(0.95, mm, nn-mm)),
      mean(containsTruth < (mm*(nn-1))/(n*(nn-mm))*qf(0.99, mm, nn-mm)),
      mean(containsTruth < (mm*(nn-1))/(n*(nn-mm))*qf(0.999, mm, nn-mm)))
