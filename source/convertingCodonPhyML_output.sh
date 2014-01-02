#!/bin/bash

die () {
  echo >&2 "$@"
  exit 1
}

## check whether input parameters are correct and input directory exists
[ "$#" -eq 2 ] || die "2 argument required: 'Phylome Data directory' & 'CDS alignment extension', $# provided"
[ -d $1 ] || die "Directory '$1' does not exist"

clear
prev=$(pwd)
printf "Changing current directory to... '%s'\n" $1;
cd $1; ext=$2;

n=0;m=0;k=0;
for i in $(ls); do
  for j in $(ls $i); do
    if [ ! -f $i/$j/$j.alg.$ext"_codonphyml_tree.txt_codon" ]; then
      n=$(expr $n + 1); continue;
    fi;

    if [ $(ls -l $i/$j/$j.alg.$ext"_codonphyml_tree.txt_codon") -eq 0 ]; then
      m=$(expr $m + 1); continue;
    fi;

    mv $i/$j/$j.alg.$ext"_codonphyml_stats.txt_codon" $i/$j/$j.tree.ml.GY_F3X4.st;
    mv $i/$j/$j.alg.$ext"_codonphyml_tree.txt_codon" $i/$j/$j.tree.ml.GY_F3X4.nw;
    lk=$(grep "^. Log-likelihood:" $i/$j/$j.tree.ml.GY_F3X4.st | awk '{ print $3}');

    printf "%s\t%s\n" "GY_F3X4" $lk > $i/$j/$j.tree.rank.ml;
    k=$(expr $k + 1);
  done;
  printf "\r%s\tFinished [%d] Unfinished [%d] NotFound [%d]" $i $k $m $n;
done;
printf "Finished [%d] Unfinished [%d] NotFound [%d]\n\n" $k $m $n;

printf "Going back to previous directory... '%s'\n" $prev;
cd $prev;
