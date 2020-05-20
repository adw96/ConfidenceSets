# README

This folder contains a number of programs that may be useful for finding the log map of a tree.

Please do write to me with a request for a program that you would like to see! I can’t help you if I don’t know what you need :)

# USAGE

** Geodesic: `java -jar geodesic.jar arg1 arg2`
     Gives the geodesic path starting at arg1 and ending at arg2

** Log map: `java -jar logmap.jar arg1 arg2`
     Gives the log map from arg1 to arg2

** Frechet mean: `java -jar meantree.jar arg1 0.0001 1000`
    Finds the sample Frechet mean tree of the trees in arg1 using the Bacak algorithm to within arg2 stepwise change running at least arg3 iterations

# EXAMPLES

`java -jar geodesic.jar tree_owen_provan.txt  tree_prime_owen_provan.txt`
`java -jar logmap.jar tree_owen_provan.txt tree_prime_owen_provan.txt`
`java -jar meantree.jar the_sample_of_trees_tmp.txt 0.0001 1000`
