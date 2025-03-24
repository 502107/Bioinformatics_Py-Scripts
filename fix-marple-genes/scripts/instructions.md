# Step1: Using new-check-primer-results.ipynb
To check the primers correspond to the correct genes, and create the reference assembly and annotations for those genes.
- Use `pass_primers-pst-2024.xlsx` that contains the primer pairs that amplified in the lab
- The code checks the files with the primer pairs designed using the Nextflow pipeline and updates the missing fields from the excel file -- output: `edited-pass_primers-pst-2024.csv`
- The merged file is then checked to make sure the gene names are correct -- output: `updated_primers_new-code.csv`
- The gff and fna files that contain the amplified and filtered genes are then created.

# Step2: Using rename_pst_metadata.ipynb
To create the primer and sample metadata files that will go into marple.
- The `marple_isolates_pst.csv` contains the accession and tree name of the new samples. The code takes the old `pst_sample_metadata.xlsx` from the marple directory and changes the old accessions to the new ones.
- The `updated_primers_new-code.csv` file, created in the previous step, is used to extract the position of each primer wrt the gene (in this case we have 1 kb padding on either side, so the code removes the overhangs and trims the sequence to just the gene length) -- output: `pst_primer_metadata.csv`

# Step3: Using search_old_genes.ipynb
To check which genes from PST130 we created primers for using Pst104e
- Output1: `old-to-new-genes.csv` has the matches from BLAST search of old genes against new reference and keeping matches with pident > 97 and % gene length match between 90-110%.
- Output2: `primers-for-old-to-new-genes.csv` has information on which genes we designed primers for.
