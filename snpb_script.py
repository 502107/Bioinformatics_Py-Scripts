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

    snp_per_base_list = []

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
        snp_per_base = snps / reference_length
        snp_per_base_list.append(snp_per_base)

    if not snp_per_base_list:
        print("No valid sequences to compare with the reference.")
        return

    # Calculate the average SNPs per base across all accessions
    average_snp_per_base = sum(snp_per_base_list) / len(snp_per_base_list)
    variance = sum((snp_per_base - average_snp_per_base) ** 2 for snp_per_base in snp_per_base_list) / len(snp_per_base_list)
    std_dev = math.sqrt(variance)

    print(f"Average SNPs per base (SNP/b) for gene '{target_gene}' using '{reference_accession}' as the reference: {average_snp_per_base:.4f}")
    print(f"STDEV SNPs per base (SNP/b) for gene '{target_gene}' using '{reference_accession}' as the reference: {std_dev:.4f}")

def calculate_snps_per_base_for_isols(cwd, target_gene, reference_accession, output_csv):
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

    if not os.path.exists(output_csv):
        with open(output_csv, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Isolate', target_gene])

    accession_data = {}

    with open(output_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row:
                accession_data[row['Isolate']] = row


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
        snp_per_base = snps / reference_length
        
        if accession in accession_data:
            accession_data[accession][target_gene] = snp_per_base
        else:
            accession_data[accession] = {'Isolate': accession, target_gene: snp_per_base}
    
    with open(output_csv, 'w', newline='') as f:
        fields = ["Isolate"] + sorted({field for row in accession_data.values() for field in row.keys() if field != "Isolate"})
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in accession_data.values():
            writer.writerow(row)

cwd = '/jic/scratch/groups/Diane-Saunders/main/MARPLE/PST/consensus'
genes_of_interest = ['jgi.p|Pucstr1|474','jgi.p|Pucstr1|9243','jgi.p|Pucstr1|20937','jgi.p|Pucstr1|11781'] # Pst425,PST130_09225,PST130_06515,PST130_14067
# genes_of_interest = ['jgi.p|Pucstr1|474','jgi.p|Pucstr1|9243']
target_gene = 'jgi.p|Pucstr1|9243'
reference_accession = '21-0005'

# for target_gene in genes_of_interest:
#     calculate_snps_per_base_for_gene(cwd, target_gene, reference_accession)
for target_gene in genes_of_interest:
    calculate_snps_per_base_for_isols(cwd, target_gene, reference_accession,'snp_per_base-per_isol.csv')



