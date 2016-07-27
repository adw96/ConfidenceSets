#!/bin/sh

#  zika_shell.sh
#  
#
#  Created by adw96 on 7/19/16.
#
# A script to populate Table 1  in my confidence set paper
# Generates 1000 sets of trees of 20, 50, 100 trees each, and assesses the coverage of a confidence set built based on this collection
# This particular script mimics the data structure of the Zika dataset

##########################################
####### START WITH 100 trees each #######

## Zika sequences had about 10,000 base pairs; but this is impractical here so we use 500

## First consider 1000 samples of 100 trees = 100,000 trees
rm -r generated100
mkdir generated100
mkdir generated100/split
mkdir generated100/split/trees

/Users/adw96/Documents/Phylogenetics/Tree\ Inference/Tree\ Inference\ Analysis/Coverage/Seq-Gen.v1.3.3/source/seq-gen -mHKY -l500 -n100000 < base_tree_zika.txt > generated100/zika_100000.dat

### work with 1000 samples of 100 trees

split -l 1800 generated100/zika_100000.dat generated100/split/zika_100000
for i in generated100/split/zika*; do split -l 9 "$i" "$i"; rm "$i"; done


### build the trees, this step takes a long time
### do these in command line
for j in generated100/split/*; do sed -i '' 's/asia_early/asia_early /g'  "$j"; /Users/adw96/Documents/Phylogenetics/PhyML/PhyML-3.1_macOS-MountainLion -i "$j"; done
cd generated100/split; for i in *_tree.txt; do mv "$i" trees/"$i"; done; cd ..; cd ..


mkdir generated100/trees

##########################################
####### 50 trees each
##########################################

# run 50 trees part!

rm -r generated50
mkdir generated50
mkdir generated50/split

/Users/adw96/Documents/Phylogenetics/Tree\ Inference/Tree\ Inference\ Analysis/Coverage/Seq-Gen.v1.3.3/source/seq-gen -mHKY -l350 -n50000 < base_tree_zika.txt > generated50/zika_50000.dat

split -l 1800 generated50/zika_50000.dat generated50/split/zika_50000
for i in generated50/split/zika*; do split -l 9 "$i" "$i"; rm "$i"; done

for j in generated50/split/*; do sed -i '' 's/asia_early/asia_early /g'  "$j"; /Users/adw96/Documents/Phylogenetics/PhyML/PhyML-3.1_macOS-MountainLion -i "$j"; done

mkdir generated50/split/trees
cd generated50/split; for i in *_tree.txt; do mv "$i" trees/"$i"; done; cd ..; cd ..

mkdir generated50/trees


##########################################
####### 20 trees each
##########################################

rm -r generated20
mkdir generated20
mkdir generated20/split

/Users/adw96/Documents/Phylogenetics/Tree\ Inference/Tree\ Inference\ Analysis/Coverage/Seq-Gen.v1.3.3/source/seq-gen -mHKY -l350 -n20000 < base_tree_zika.txt > generated20/zika_20000.dat

split -l 1800 generated20/zika_20000.dat generated20/split/zika_20000
for i in generated20/split/zika*; do split -l 9 "$i" "$i"; rm "$i"; done

for j in generated20/split/*; do sed -i '' 's/asia_early/asia_early /g'  "$j"; /Users/adw96/Documents/Phylogenetics/PhyML/PhyML-3.1_macOS-MountainLion -i "$j"; done

mkdir generated20/split/trees
cd generated20/split; for i in *_tree.txt; do mv "$i" trees/"$i"; done; cd ..; cd ..
mkdir generated20/trees

#sed -i '' 's/asia_early/asia_early /g' generated20/split/zika_20000aaak
#
#sed -e "/asia_early/a\\
#''" < generated20/split/zika_20000aaae
#
#sed '/asia_early/a\asia\' generated20/split/zika_20000aaae
#
#perl -pi -i -e '$_ .= qq(asia_early\n) if /asia_early/' generated20/split/zika_20000aaae
#
#sed 's/asia_early/&\
#/g' generated20/split/zika_20000aaaf
#
#sed 's/asia_early/\
#&/g' generated20/split/zika_20000aaag
#
#perl -lne 'print $_;print "Cygwin" if(/Fedora/);' generated20/split/zika_20000aaah
#perl -lne 'print $_;print "\n" if(/asia_early/);' generated20/split/zika_20000aaah
#
#sed -i '' 's/asia_early/asia_early/\n/g' generated20/split/zika_20000aaaj
#
#sed -i -e '/asia_early/a\n' generated20/split/zika_20000aaaj
#
#sed "s/regexp/\\`echo -e '\n\r'`/g"
#
#sed -i '' 's/asia_early/asia_early n/g' generated20/split/zika_20000aaaj
#
#sed -i '' 's/asia_early/asia_early /g' generated20/split/zika_20000aaak
#

#sed -i '' 's/asia_early/asia_early /g' generated20/split/zika_20000aaaa


#/Users/adw96/Documents/Phylogenetics/PhyML/PhyML-3.1_macOS-MountainLion -i generated20/split/zika_20000aaaz

#/Users/adw96/Documents/Phylogenetics/PhyML/PhyML-3.1_macOS-MountainLion -i generated20/split/zika_20000aaaa

