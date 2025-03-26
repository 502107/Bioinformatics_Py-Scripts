#!/bin/bash
#SBATCH -J IQtree
#SBATCH --mem 100G
#SBATCH --nodes=2                     # Number of nodes
#SBATCH --ntasks-per-node=32          # Number of tasks (processors)
#SBATCH -t 02:00:00
#SBATCH -p jic-long,nbi-long,RG-Diane-Saunders,jic-short
#SBATCH -o IQTree.%j.out
#SBATCH -e IQTree.%j.err

source iqtree-2.0.5

while getopts ":i:h" opt; do
  case $opt in
    i) input="$OPTARG"
    ;;
    h) echo "Usage: sbatch $0 -i <input file (.phy)> "
       exit 0
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
       exit 1
    ;;
  esac
done

if [ -z "$input" ]; then
    echo "Input file must be provided. Use -h for help."
    exit 1
fi

iqtree2 -s "$input" -st DNA -m GTR+I+G -T 2 -B 1000 


<<comment
corresponding raxml-ng script:
raxml-ng --all --msa PST_concatinatelan_serbia.phy --model GTR+G --tree pars{10},rand{10} --prefix pst_serbia --bs-trees 100 --seed 2  

parameters used here

-T 2 -- Use 2 CPUs
-m GTR+I+G -- Infer maximum likelihood tree using GTR+I+G
-B 1000 -- Replicates for ultrafast bootstrap (>=1000)

comment

