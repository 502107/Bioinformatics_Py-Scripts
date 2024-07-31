from Bio import SeqIO
import re
import os

REF_dir = '/jic/scratch/groups/Diane-Saunders/loizos/PGT/REF/CRL75/GCA_000149925.1/'
ref_ann = os.path.join(REF_dir,'genomic.gff')
ref_fna = os.path.join(REF_dir,'GCA_000149925.1_ASM14992v1_genomic.fna')
gene_list = '276genes.txt'

gene_dict = {}

with open(gene_list, 'r') as f:
    genes = f.read().splitlines()

with open(ref_ann, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            ann = line.split('\t')
            if ann[2] == 'gene':
                gene_id = re.search('ID=gene-(.*?);', ann[8]).group(1)
                if gene_id in genes:
                    scaff = ann[0]
                    start = int(ann[3])
                    end = int(ann[4])
                    gene_dict[gene_id] = [scaff, start, end]
                    

with open(ref_fna, 'r') as f, open('genes.fasta', 'w') as outf:
    for record in SeqIO.parse(f, 'fasta'):
        for gene_id, [scaff, start, end] in gene_dict.items():
            if record.id == scaff:
                gene_seq = record.seq[start:end]
                print('>' + gene_id, file=outf)
                print(gene_seq, file=outf)

