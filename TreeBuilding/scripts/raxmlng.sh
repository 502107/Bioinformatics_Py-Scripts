#!/bin/bash
#SBATCH -J raxml_tree
#SBATCH --mem 100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -t 02:00:00
#SBATCH -p jic-long,nbi-long,RG-Diane-Saunders,jic-short
#SBATCH -o pst_serbia_tree.%j.out
#SBATCH -e pst_serbia_tree.%j.err

source package ource package /tgac/software/production/bin/raxml-ng-0.8.0

while getopts ":i:o:h" opt; do
  case $opt in
    i) input="$OPTARG"
    ;;
    o) output_pref="$OPTARG"
    ;;
    h) echo "Usage: sbatch raxml.sh -i <input file> -o <output prefix>"
       exit 0
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
       exit 1
    ;;
  esac
done

if [ -z "$input" ] || [ -z "$output_pref" ]; then
    echo "Both input file and output prefix must be provided. Use -h for help."
    exit 1
fi

raxml-ng --all --msa "$input" --model GTR+G --tree pars{10},rand{10} --prefix "$output_pref" --bs-trees 100 --seed 2 

##raxmlHPC-PTHREADS-SSE3: run RAxML using the PTHREADS-SSE3 version, which supports parallel processing on multiple threads.

#-T 20: number of threads to be used for parallel processing, set to 20 threads.

#`-s` option specifies the input alignment file in PHYLIP format

#`-m` option is used to specify the substitution model to be used during the analysis. In this case, "GTRGAMMA" represents the General Time Reversible (GTR) model with a gamma distribution of rates across sites.

#`-n` option specifies the name of the output file that will contain the resulting tree(s) in Newick format

# `-p` option is used to specify the random number seed for reproducibility

#-N 20`: This option specifies the number of starting trees to be generated during the analysis. By default, RAxML generates a single starting tree

