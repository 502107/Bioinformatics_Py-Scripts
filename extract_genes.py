import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def parse_gff_genes(gff_file):
    """Parse GFF file to extract gene coordinates and attributes"""
    genes = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            feature_type = fields[2]
            if feature_type != 'gene':
                continue
                
            contig = fields[0]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            gene_id = None
            gene_name = None
            
            for attr in attributes.split(';'):
                if attr.startswith('ID='):
                    gene_id = attr.split('=')[1]
                elif attr.startswith('Name='):
                    gene_name = attr.split('=')[1]
                elif attr.startswith('gene_id='):
                    gene_id = attr.split('=')[1]
                elif attr.startswith('gene_name='):
                    gene_name = attr.split('=')[1]
            
            if gene_id:
                genes[gene_id] = {
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'name': gene_name or gene_id
                }
    
    return genes

def extract_gene_sequences(fasta_file, genes, output_file):
    """Extract gene sequences from fasta file"""
    ref_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    
    gene_records = []
    
    for gene_id, gene_info in genes.items():
        contig = gene_info['contig']
        start = gene_info['start']
        end = gene_info['end']
        strand = gene_info['strand']
        name = gene_info['name']
        
        if contig not in ref_sequences:
            print(f"Warning: Contig {contig} not found in reference fasta")
            continue
        
        seq = ref_sequences[contig].seq[start:end]
        
        if strand == '-':
            seq = seq.reverse_complement()
        
        record = SeqRecord(
            seq,
            id=name,
            description=f"gene_id={gene_id} {contig}:{start+1}-{end} strand={strand}"
        )
        
        gene_records.append(record)
    
    SeqIO.write(gene_records, output_file, 'fasta')
    print(f"Extracted {len(gene_records)} genes to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Extract gene sequences from reference assembly')
    parser.add_argument('--fasta', required=True, help='Reference assembly fasta file')
    parser.add_argument('--gff', required=True, help='GFF annotation file')
    parser.add_argument('--output', required=True, help='Output fasta file with gene sequences')
    
    args = parser.parse_args()
    
    genes = parse_gff_genes(args.gff)
    
    extract_gene_sequences(args.fasta, genes, args.output)

if __name__ == '__main__':
    main()
