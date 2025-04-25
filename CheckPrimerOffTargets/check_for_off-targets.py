# This script checks the diversity of primers against a modified reference genome.
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
import os

ref_fna = 'assets/Triticum_aestivum.IWGSC.dna.toplevel.fa'
ref_gff = 'assets/Triticum_aestivum.IWGSC.60.gff3'

primer_seqs_F = ['GCCACGACCACTTCTCTAGC', 'TTCTCTAGCTCATGCCCGTC', 'CGCCTGTGTCCCGATAGATC', 
                 'AAGACATGGACGAGCAGTGG', 'AAGAATCACGAGGCTGAGGC', 'GCAGGGTCTCCAAGCTCTAC', 
                 'CTCCAAGCTCTACGTCCACG', 'TACCAATGGCGGAAGTACGG']
primer_seqs_R = ['AAGTCTGCCCGATCATCCAC', 'GGAAGTCTGCCCGATCATCC', 'ACATGAAGTCCTCCTCCACG', 
                 'CTGAGCCTCCTGTTCTCCTC', 'ACGACTGGTGGTGGTTGTTG', 'TGCACCTTCTTCTTCACCGG', 
                 'CGCTGCACCTTCTTCTTCAC', 'CACCAGAGCCTTGTTGCTTG']
target_gene = 'TraesCS7D02G497300.1'

def search_mismatch(primer, sequence, max_mismatch=4, forward=True):
    """
    Return True if the primer is diverse enough (at least 1 SNP in the last 3bp).
    """
    primer_len, seq_len = len(primer), len(sequence)
    assert seq_len >= primer_len, 'Primer is larger than the sequence'

    if primer in sequence:  # Exact match - bad primer
        # Print +/- 50 bp around the sequence match
        location = sequence.find(primer)
        start = max(0, location - 50)
        end = min(seq_len, location + len(primer) + 50)
        return False, sequence[start:end]

    for i in range(seq_len - primer_len + 1):
        segment = sequence[i:i + primer_len]
        mm = sum(1 for a, b in zip(primer, segment) if a != b)
        if mm < max_mismatch:
            if forward:
                last_mm = sum(1 for a, b in zip(primer[-3:], segment[-3:]) if a != b)
            else:
                last_mm = sum(1 for a, b in zip(primer[:3], segment[:3]) if a != b)
            if last_mm > 0:
                return True, None
            else:
                # Print +/- 50 bp around the sequence match
                location = sequence.find(segment)
                start = max(0, location - 50)
                end = min(seq_len, location + len(segment) + 50)
                return False, segment

    return True, None


def offtarget_in_host(primer, host_genome, forward=True):
    """
    Check for off-targets in the host genome.
    """
    host_genome = str(Seq(host_genome))
    primer_rc = str(Seq(primer).reverse_complement())

    if forward:
        is_diverse, segment = search_mismatch(primer, host_genome)
    else:
        is_diverse, segment = search_mismatch(primer_rc, host_genome, forward=False)

    return not is_diverse, segment

def modify_sequence(ref_fna, ref_gff, target_gene):
    """
    Modify the reference genome by masking the target gene with 'N'.
    """
    if not os.path.exists(ref_fna):
        print(f"Error: Reference genome file '{ref_fna}' not found.")
        return None

    if not os.path.exists(ref_gff):
        print(f"Error: GFF file '{ref_gff}' not found.")
        return None

    modified_sequence = ""
    target_gene_loc = {}
    found = False

    with open(ref_gff, "r") as gff_file:
        for line in gff_file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "mRNA" and "ID=" in fields[8]:
                gene_id = fields[8].split(':')[1].split(";")[0]
                if gene_id == target_gene:
                    chrom = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    target_gene_loc[gene_id] = [(chrom, start, end)]
                    found = True

    if not found:
        print(f"Target gene {target_gene} not found in GFF file.")
        return None

    for record in SeqIO.parse(ref_fna, "fasta"):
        for gene_id, location in target_gene_loc.items():
            chrom, start, end = location[0]
            if record.id == chrom:
                sequence = str(record.seq)
                sequence = sequence[:start - 1] + "N" * (end - start + 1) + sequence[end:]
                modified_sequence += sequence
            else:
                modified_sequence += str(record.seq)

    return modified_sequence

def main():
    modified_sequence = modify_sequence(ref_fna, ref_gff, target_gene)
    if modified_sequence is None:
        return

    diverse_primers = []
    non_diverse_primers = []

    for primer in tqdm(primer_seqs_F, desc="Checking forward primers"):
        has_offtarget, segment = offtarget_in_host(primer, modified_sequence, forward=True)
        if has_offtarget:
            print(f"The forward primer {primer} is not diverse enough.")
            if segment:
                print(f"Location: {segment}")
            non_diverse_primers.append(primer)
        else:
            print(f"The forward primer {primer} is diverse enough.")
            diverse_primers.append(primer)

    for primer in tqdm(primer_seqs_R, desc="Checking reverse primers"):
        has_offtarget, segment = offtarget_in_host(primer, modified_sequence, forward=False)
        if has_offtarget:
            print(f"The reverse primer {primer} is not diverse enough.")
            if segment:
                print(f"Location: {segment}")
            non_diverse_primers.append(primer)
        else:
            print(f"The reverse primer {primer} is diverse enough.")
            diverse_primers.append(primer)

    print("\nSummary:")
    print(f"Diverse primers: {len(diverse_primers)}")
    print(f"Non-diverse primers: {len(non_diverse_primers)}")

if __name__ == "__main__":
    main()
