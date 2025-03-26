#!/bin/bash
#SBATCH -J fasttree
#SBATCH --mem 100G
#SBATCH --nodes=2                     # Number of nodes
#SBATCH --ntasks-per-node=32          # Number of tasks (processors)
#SBATCH -t 10:00:00
#SBATCH -p jic-long,nbi-long
#SBATCH -o FastTree.%j.out
#SBATCH -e FastTree.%j.err

source fasttree-2.1.9

while getopts ":i:o:h" opt; do
  case $opt in
    i) input="$OPTARG"
    ;;
    o) output="$OPTARG"
    ;;
    h) echo "Usage: sbatch fasttree.sh -i <input file (.fa)> -o <output file (.nwk)>"
       exit 0
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
       exit 1
    ;;
  esac
done

if [ -z "$input" ] || [ -z "$output" ]; then
    echo "Both input and output files must be provided. Use -h for help."
    exit 1
fi

fasttree -gtr -nt < "$input" > "$output"


<<comment
corresponding raxml-ng script:
raxml-ng --all --msa PST_concatinatelan_serbia.phy --model GTR+G --tree pars{10},rand{10} --prefix pst_serbia --bs-trees 100 --seed 2  

parameters used here

-gtr -- generalized time-reversible model (nucleotide alignments only)
-nt -- working with nucleotide sequences.
-mlacc 2 -slownni -- make the maximum-likelihood NNIs more exhaustive (~4x slower for 5,000 proteins)

FYI:

FastTree does not have built-in support for specifying bootstrap trees in the same way that RAxML does. 
In FastTree, the bootstrap analysis is typically performed internally, and it doesnt provide options for specifying the number of bootstrap replicates or controlling how the bootstrap trees are generated.

FastTree is designed to be a fast and efficient tool for estimating phylogenetic trees, particularly for large datasets.
It uses a different approach than RAxML for branch support calculation. FastTree employs local support values rather than traditional bootstrap values.
These local support values are based on the Shimodaira-Hasegawa test and can provide an assessment of the statistical support for each branch in the tree.

[REF: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2835736/]

comment

