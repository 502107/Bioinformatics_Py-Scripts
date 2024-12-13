import os
import re
import sys
import getopt
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO

def nt_to_aa(input):
    fna=input
    aa_out=f'{fna.split(".")[0]}_aa.fna'
    for record in SeqIO.parse(fna,"fasta"):
        nt_seq = Seq(record.seq)
        aa_seq = nt_seq.translate()
        with open(aa_out,'a') as fout:
            fout.write(f'>{record.id}\n{aa_seq}\n')
            
def blast(query, subject, gff, run='megablast'):
    subprocess.run(["blastn", "-query", query, "-subject", subject, "-out", "blast.out", "-outfmt", "6 qseqid sseqid pident bitscore sstart send sseq", "-task", run], check=True)
    
    if gff:
        sseqid_dict = {}
        locus_tag_pattern = r'locus_tag=([^;]*)'
        id_pattern = r'ID=([^;]*)'
    
        with open('blast.out', 'r') as fb:
            for bline in fb:
                bline = bline.strip().split('\t')
                sseqid = bline[1]
                sstart = int(bline[4])
                send = int(bline[5])
                with open(gff, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        line = line.strip().split('\t')
                        if line[2] == 'CDS' and line[0] == sseqid:
                            gene_loc = line[3:5]
                            if sseqid == line[0]:
                                if (sstart >= int(gene_loc[0]) and send <= int(gene_loc[1])) or (sstart <= int(gene_loc[0]) and send >= int(gene_loc[1])):
                                    # Get gene name by locus_tag or ID in gff
                                    locus_tag_match = re.search(locus_tag_pattern, line[8])
                                    id_match = re.search(id_pattern, line[8])
                                    if locus_tag_match:
                                        gene_name = locus_tag_match.group(1)
                                    elif id_match:
                                        gene_name = id_match.group(1)
                                    else:
                                        gene_name = line[8]
                                    if bline[0] not in sseqid_dict:
                                        sseqid_dict[bline[0]] = gene_name

        for sseqid, gene in sseqid_dict.items():
            print(f'Potential gene {gene} found for BLAST match in {sseqid}')

def get_gene(gff, fna, gene):
    gene_list = {}
    gene_count = {}
    gene_seqs = {}
    og_gene = gene

    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if line[2] == 'CDS':
                gene_loc = line[3:5]
                sseqid = line[0]
                sstart, send = int(gene_loc[0]), int(gene_loc[1])
                gene_name = line[8]
                try:
                    if gene in gene_name:
                        gene_count[gene] = gene_count.get(gene,0) + 1
                        if gene not in gene_list:
                            gene_list[gene] = []
                        gene_list[gene].append((sseqid, sstart, send))
                except:
                    print(f'Gene {gene} not found in GFF file')

    for gene,locs in gene_list.items():
        # Sort by sstart
        locs.sort(key=lambda x: x[1])
        for i, loc in enumerate(locs):
            if loc[1] <= loc[2]:
                gene_suffix = chr(ord('A') + i)
                gene_seqs[gene + gene_suffix] = loc
            else:
                gene_suffix = chr(ord('A') + i)
                gene_seqs[gene + gene_suffix] = loc[0], loc[2], loc[1]

    combined_seq = ''
    for gene, gene_loc in gene_seqs.items():
        sseqid = gene_loc[0]
        sstart = gene_loc[1]
        send = gene_loc[2]
        for record in SeqIO.parse(fna,"fasta"):
            if sseqid == record.id:
                combined_seq += str(record.seq[sstart-1:send])
    
    with open(f'{og_gene}.fna','w') as fout:
        fout.write(f'>{og_gene}\n{combined_seq}\n')
            
def align(inputs, run=False):
    if os.path.exists('alignment.fna'):
        os.remove('alignment.fna')
    for input in inputs:
        for record in SeqIO.parse(input,"fasta"):
            with open('alignment.fna','a') as fout:
                fout.write(f'>{record.id}_{input.split(".")[0]}\n{record.seq}\n')
    if run == 'nt':
        subprocess.run(["clustalw2", "-infile=" + "alignment.fna", "-outfile=" + 'alignment.aln', "-output=fasta"], check=True)
    elif run == 'aa':
        subprocess.run(["clustalw2", "-infile=" + "alignment.fna", "-outfile=" + 'alignment.aln', "-output=fasta", "-type=Protein"], check=True)

def extract_cds(gff, fna):
    gene_list = {}
    gene_count = {}
    gene_seqs = {}

    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if line[2] == 'CDS':
                gene_loc = line[3:5]
                sseqid = line[0]
                sstart, send = int(gene_loc[0]), int(gene_loc[1])

                attributes = line[8].split(';')
                gene_name = [attr.split('=')[1] for attr in attributes if attr.startswith('ID')][0]

                gene_count[gene_name] = gene_count.get(gene_name, 0) + 1

                if gene_count[gene_name] == 1:
                    gene_seqs[gene_name] = ''

    with open(fna) as f:
        for line in f:
            if line.startswith('>'):
                sseqid = line[1:].strip()
            else:
                if sseqid in gene_seqs:
                    gene_seqs[sseqid] += line.strip()

    with open(f"{gff.split('.')[0]}_CDS.fasta", 'w') as f:
        for gene_name, sequence in gene_seqs.items():
            f.write(f'>{gene_name}\n{sequence}\n')


def main(argv):
    task = ''
    input = ''
    gff = False
    fna = ''
    gene = ''
    query = ''
    subject = ''
    run = False
    a_inputs = []
    try:
        opts, args = getopt.getopt(argv,"ht:i:g:f:c:1:2:3:4:r:q:s:",["task=","input=","gff=","fna=","gene=","input1=","input2=","input3=","input4=","run=","query=","subject="])
    except getopt.GetoptError:
        print(f'Usage: python bioinformatic_tools -t <task>\nTasks:\npython bioinformatic_tools -t nt_to_aa -i <input_fna>\npython bioinformatic_tools -t get_gene -g <input_gff> -f <input_fna> -c <gene_name>\npython bioinformatic_tools -t align -1 <input_fna> -2 <input_fna> [up to 4 input sequences; optional: -r [nt/aa]]\npython bioinformatic_tools -t blast -q <query> -s <subject> [optional: -r [blastn/megablast] -g <annotation_gff>]\npython bioinformatic_tools -t extract_cds -g <input_gff> -f <input_fna>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(f'Usage: python bioinformatic_tools -t <task>\nTasks:\npython bioinformatic_tools -t nt_to_aa -i <input_fna>\npython bioinformatic_tools -t get_gene -g <input_gff> -f <input_fna> -c <gene_name>\npython bioinformatic_tools -t align -1 <input_fna> -2 <input_fna> [up to 4 input sequences; optional: -r [nt/aa]]\npython bioinformatic_tools -t blast -q <query> -s <subject> [optional: -r [blastn/megablast] -g <annotation_gff>]\npython bioinformatic_tools -t extract_cds -g <input_gff> -f <input_fna>')
            sys.exit()
        elif opt in ("-t", "--task"):
            task = arg
        elif opt in ("-i", "--input"):
            input = arg
        elif opt in ("-g", "--gff"):
            gff = arg
        elif opt in ("-f", "--fna"):
            fna = arg
        elif opt in ("-c", "--gene"):
            gene = arg
        elif opt in ("-1", "--input1"):
            a_inputs.append(arg)
        elif opt in ("-2", "--input2"):
            a_inputs.append(arg)
        elif opt in ("-3", "--input3"):
            a_inputs.append(arg)
        elif opt in ("-4", "--input4"):
            a_inputs.append(arg)
        elif opt in ("-r", "--run"):
            run = arg
        elif opt in ("-q", "--query"):
            query = arg
        elif opt in ("-s", "--subject"):
            subject = arg
    if task == 'nt_to_aa':
        nt_to_aa(input)
    elif task == 'extract_cds':
        extract_cds(gff, fna)
    elif task == 'get_gene' and gff != False:
        get_gene(gff,fna,gene)
    elif task == 'align':
        if run == 'nt' or run == 'aa':
            align(a_inputs, run=run)
        else:
            align(a_inputs)
    elif task == 'blast':
        if run == 'blastn':
            blast(query, subject, gff, run='blastn')
        elif run == 'megablast' or run == False:
            blast(query, subject, gff)

if __name__ == "__main__":
   main(sys.argv[1:])
   
