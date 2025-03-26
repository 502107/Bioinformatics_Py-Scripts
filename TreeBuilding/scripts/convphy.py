# last update: 22/11/2023

import sys
import os

# User input (filename) processing

if len(sys.argv) != 2:
    print("Provide filename of a concatinated phylip file: python convphy.py <.phy>")
    sys.exit(1)

input_file=sys.argv[1]
filename = os.path.basename(input_file)

inp_noext=filename.split('.')[0]


with open(input_file, 'r') as pg, open(inp_noext+'.fa','w') as cpg, open(inp_noext+'_onlySNP.fa','w') as snpg, open(inp_noext+'_onlySNP.phy','w') as snphy:

    next(pg)

    head_ls = []
    seq_ls = []
    head_len=0

    for line in pg:
        chars = line.split()
        head = chars[0]
        seq = chars[1]

        head_ls.append(head)
        seq_ls.append(seq)

        cpg.write('>' + head + '\n' + seq + '\n')

	head_len +=1

    positions = ''
    for pos in range(len(seq_ls[0])):

        # Get the nucleotides at each position
        nuc_set = set(seq_ls[i][pos] for i in range(len(seq_ls)) if seq_ls[i][pos] in 'ACGTU')

        # Check if the set has more than one element, ie. there is a SNP
        if len(nuc_set) > 1:
            positions += str(pos + 1) + ' '

    for i in range(len(head_ls)):
        snpg.write('>' + head_ls[i] + '\n')
        snphy.write(head_ls[i]+ ' ')
	seq_len=0
        for pos in positions.split():
            # Write the nucleotide at each position to the SNP output file
            snpg.write(seq_ls[i][int(pos) - 1])
            snphy.write(seq_ls[i][int(pos) - 1])
	    seq_len +=1
        snpg.write('\n')
        snphy.write('\n')

with open(inp_noext+'_onlySNP.phy', 'r') as rtmp, open('tmp', 'w') as wtmp:
    wtmp.write(str(head_len) + ' ' + str(seq_len)+ '\n' + rtmp.read())
    os.remove(inp_noext+'_onlySNP.phy')
    os.rename('tmp',inp_noext+'_onlySNP.phy')

