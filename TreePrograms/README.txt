README

This folder contains a number of programs that may be useful for finding the log map of a tree. 

Please do write to me with a request for a program that you would like to see! I can’t help you if I don’t know what you need :)

USAGE

** Geodesic: java -jar geodesic.jar arg1 arg2
     Gives the geodesic path starting at arg1 and ending at arg2

** Log map: java -jar logmap.jar arg1 arg2
     Gives the log map from arg1 to arg2

** Log map of a dataset: java -jar map_dataset.jar arg1 (arg2) (arg3) (arg4)
     Calculates the Frétchet mean of the trees in arg1 to within a threshold of arg2 or arg3 iterations (whichever happens first), and writes the mean tree to output or to arg4. 
     I like to write my trees into R, so I often use “java -jar map_dataset.jar test1.txt > my_script.R”

EXAMPLES

java -jar geodesic.jar tree_owen_provan.txt  tree_prime_owen_provan.txt 
java -jar logmap.jar tree_owen_provan.txt tree_prime_owen_provan.txt 
java -jar map_dataset.jar test1.txt 100 10 out.txt
java -jar map_dataset.jar test1.txt 100 10 
java -jar map_dataset.jar test1.txt 