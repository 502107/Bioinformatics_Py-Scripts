import os
import csv
import math

def calculate_snps_per_base_for_gene(cwd, target_gene, reference_accession):
    seq_dict = {}

    for file in os.listdir(cwd):
        file_path = os.path.join(cwd, file)
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if target_gene in line:
                        isol = file.split('_')[-4]  # Extract isolate ID from filename
                        seq_dict[isol] = next(f).strip()  # Store the sequence
                        break
        except OSError as e:
            print(f"Error opening file {file_path}: {e}")

    if reference_accession not in seq_dict:
        print(f"Error: Reference accession '{reference_accession}' not found in the sequences.")
        return

    # Use a specified accession's sequence as the reference
    reference = seq_dict[reference_accession]
    reference_length = len(reference)

    snp_counts = []
    sequence_count = 0

    for accession, seq in seq_dict.items():
        if accession == reference_accession:
            continue  # Skip the reference sequence itself
        if len(seq) != reference_length:
            print(f"Warning: Sequence length mismatch for accession '{accession}'. Skipping sequence.")
            continue
        amb_count = seq.count('?')
        amb_pct = amb_count / reference_length
        if amb_pct > 0.4:
            # print(f"Warning: Sequence '{accession}' has >40% ambiguities. Skipping sequence.")
            continue
        snps = sum(1 for a, b in zip(reference, seq) if a != b)
        snp_counts.append(snps)
        sequence_count += 1

    if sequence_count == 0:
        print("No valid sequences to compare with the reference.")
        return

    # Calculate the average SNPs per base
    average_snp_per_base = sum(snp_counts) / (reference_length * sequence_count)
    variance = sum((snps - average_snp_per_base * reference_length) ** 2 for snps in snp_counts) / sequence_count
    std_dev = math.sqrt(variance)

    print(f"Average SNPs per base (SNP/b) for gene '{target_gene}' using '{reference_accession}' as the reference: {average_snp_per_base:.4f}")
    print(f"STDEV SNPs per base (SNP/b) for gene '{target_gene}' using '{reference_accession}' as the reference: {std_dev:.4f}")
    
cwd = '/jic/scratch/groups/Diane-Saunders/main/MARPLE/PST/consensus'
genes_of_interest = ['jgi.p|Pucstr1|474','jgi.p|Pucstr1|9243','jgi.p|Pucstr1|20937','jgi.p|Pucstr1|11781'] # Pst425,PST130_09225,PST130_06515,PST130_14067
target_gene = 'jgi.p|Pucstr1|9243'
reference_accession = '21-0005'

for target_gene in genes_of_interest:
    calculate_snps_per_base_for_gene(cwd, target_gene, reference_accession)