#!/bin/bash

source package aeee87c4-1923-4732-aca2-f2aff23580cc
source package 6394519c-541f-479c-b064-dd0b912eac04

# Pst425+2kb flanking each side
coord="NC_063032.1:2672372-2676809"
bam=$1
acc=$(basename "$bam" .bam)
reference='GCF_021901695.1_Pst134E36_v1_pri_genomic.fna'

tmp_dir=$(mktemp -d)
trap 'rm -rf $tmp_dir' EXIT

# Step 1: Extract reads overlapping the region of interest
echo "Extracting reads from region: $coord"
samtools view -b "$bam" "$coord" > "$tmp_dir/Pst425_${acc}.bam"

# Step 2: Generate BED coverage information for the region
echo "Generating BED coverage file"
bedtools genomecov -ibam "$tmp_dir/Pst425_${acc}.bam" -bga > "$tmp_dir/coverage.bed"

# Step 3: Filter the BED file for the exact region of interest
echo "Filtering BED file for target region"
echo -e "NC_063032.1\t2672372\t2676809" > "$tmp_dir/region.bed"
bedtools intersect -a "$tmp_dir/coverage.bed" -b "$tmp_dir/region.bed" | awk '$4 > 0' > "$tmp_dir/filtered.bed"

# Step 4: Merge overlapping or adjacent intervals in the filtered BED file
echo "Merging overlapping intervals"
bedtools merge -i "$tmp_dir/filtered.bed" > "$tmp_dir/merged.bed"

# Step 5: Use the merged BED file to extract only the covered sequence
echo "Extracting sequence for the target region"
bedtools getfasta -fi "$reference" -bed "$tmp_dir/merged.bed" -fo "$tmp_dir/extracted_sequence.fna"

# Combine all fragmented sequences into a single FASTA entry
echo "Combining fragmented sequences"
awk 'BEGIN {printf(">Pst425_%s\n", "'${acc}'")} !/^>/ {printf("%s", $0)} /^>/ && NR > 1 {printf("\n")}' "$tmp_dir/extracted_sequence.fna" > "Pst425_${acc}.fna"

echo "Sequence saved to Pst425_${acc}.fna"
