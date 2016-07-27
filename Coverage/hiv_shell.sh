#!/bin/sh

#  hiv_shell.sh
#  
#
#  Created by adw96 on 7/19/16.
#
# A script to populate Table 1  in my confidence set paper
# Generates 1000 sets of trees of 20, 50, 100 trees each, and assesses the coverage of a confidence set built based on this collection
# This particular script mimics the data structure of the HIV dataset

##########################################
####### START WITH 100 trees each #######

## HIV sequences had about 350 base pairs

## First consider 1000 samples of 100 trees = 100,000 trees
rm -r generated100
mkdir generated100
mkdir generated100/split

/Users/adw96/Documents/Phylogenetics/Tree\ Inference/Tree\ Inference\ Analysis/Coverage/Seq-Gen.v1.3.3/source/seq-gen -mHKY -l350 -n100000 < hiv_base_tree.txt > generated100/hiv_100000.dat

### work with 1000 samples of 100 trees

split -l 1800 generated100/hiv_100000.dat generated100/split/hiv_100000
for i in generated100/split/hiv*; do split -l 6 "$i" "$i"; rm "$i"; done


### build the trees, this step takes a long time
### do these in command line
for i in generated100/split/*; do /Users/adw96/Documents/Phylogenetics/PhyML/PhyML-3.1_macOS-MountainLion -i "$i"; done
mkdir generated100/split/trees
cd generated100/split; for i in *_tree.txt; do mv "$i" trees/"$i"; done; cd ..; cd ..

##########################################
####### 50 trees each
##########################################

rm -r generated50
mkdir generated50
mkdir generated50/split

/Users/adw96/Documents/Phylogenetics/Tree\ Inference/Tree\ Inference\ Analysis/Coverage/Seq-Gen.v1.3.3/source/seq-gen -mHKY -l350 -n50000 < hiv_base_tree.txt > generated50/hiv_50000.dat

split -l 1800 generated50/hiv_50000.dat generated50/split/hiv_50000
for i in generated50/split/hiv*; do split -l 6 "$i" "$i"; rm "$i"; done

for i in generated50/split/*; do /Users/adw96/Documents/Phylogenetics/PhyML/PhyML-3.1_macOS-MountainLion -i "$i"; done
mkdir generated50/split/trees
cd generated50/split; for i in *_tree.txt; do mv "$i" trees/"$i"; done; cd ..; cd ..


##########################################
####### 20 trees each
##########################################

rm -r generated20
mkdir generated20
mkdir generated20/split

/Users/adw96/Documents/Phylogenetics/Tree\ Inference/Tree\ Inference\ Analysis/Coverage/Seq-Gen.v1.3.3/source/seq-gen -mHKY -l350 -n20000 < hiv_base_tree.txt > generated20/hiv_20000.dat

split -l 1800 generated20/hiv_20000.dat generated20/split/hiv_20000
for i in generated20/split/hiv*; do split -l 6 "$i" "$i"; rm "$i"; done

for i in generated20/split/*; do /Users/adw96/Documents/Phylogenetics/PhyML/PhyML-3.1_macOS-MountainLion -i "$i"; done
mkdir generated20/split/trees
cd generated20/split; for i in *_tree.txt; do mv "$i" trees/"$i"; done; cd ..; cd ..
