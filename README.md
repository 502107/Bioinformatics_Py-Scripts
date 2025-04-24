# Bioinformatics Python Scripts

This repository contains a collection of Python scripts and tools for various bioinformatics analyses. The tools are designed to help with tasks such as sequence manipulation, alignment analysis, gene extraction, and visualization.

## Dependencies

The scripts require the following dependencies:
- Python 3.x
- Biopython
- BLAST+
- ClustalW2
- SAMtools (for BAM file processing)

To install Python dependencies:
```bash
pip install biopython
```

For system tools (on macOS):
```bash
brew install blast
brew install clustalw
brew install samtools
```

## Repository Structure

- `bioinformatic_tools.py`: Main utility script with various bioinformatics functions
- `create_indices.sh`: Script for creating genome indices
- `extract_sequence_from_bam.sh`: Script for extracting sequences from BAM files
- `snpb_script.py`: Script for SNP analysis
- `Fix-alignments.ipynb`: Jupyter notebook for fixing sequence alignments

### Specialized Directories
- `Circosplt-Map_contig_to_chrom/`: Tools for mapping contigs to chromosomes and circular genome visualisation
- `FindHaplotypesIsoforms/`: Scripts for identifying haplotypes and isoforms
- `PrimerMSAplot/`: Multiple sequence alignment plotting
- `SIZ1_domains/`: Analysis tools for SIZ1 protein domains and visualisation
- `TreeBuilding/`: Scripts for phylogenetic tree construction
- `UpsetPlot/`: Tools for creating UpSet plots
- `fix-marple-genes/`: Scripts for fixing gene annotations
- `search_for_old_genes/`: Tools for finding genes from one genome assembly in another

## Main Tools and Usage

### bioinformatic_tools.py

1. **Nucleotide to Amino Acid Translation**
```bash
python bioinformatic_tools.py -t nt_to_aa -i <input_fna>
```

2. **Gene Extraction**
```bash
python bioinformatic_tools.py -t get_gene -g <input_gff> -f <input_fna> -c <gene_name>
```

3. **Sequence Alignment**
```bash
python bioinformatic_tools.py -t align -1 <input_fna1> -2 <input_fna2> [-r nt/aa]
```

4. **BLAST Analysis**
```bash
python bioinformatic_tools.py -t blast -q <query> -s <subject> [-r blastn/megablast] [-g <annotation_gff>]
```

5. **CDS Extraction**
```bash
python bioinformatic_tools.py -t extract_cds -g <input_gff> -f <input_fna>
```

### Shell Scripts

1. **Creating Genome Indices**
```bash
bash create_indices.sh --refPath <reference_genome_directory> --refGenomeName <reference_genome_fasta> --refDescriptionName <reference_gff>
```

2. **Extracting Sequences from BAM Files**
```bash
bash extract_sequence_from_bam.sh <input.bam>
```

## Specialized Directory Usage

### Circosplt-Map_contig_to_chrom
This directory contains tools for mapping contigs to chromosomes and creating circular genome visualizations.

1. **Gene Extraction and Mapping**
Extract the sequence of the genes of interest using `extract-genes.py`, from the old reference assembly.
```python
python extract-genes.py
```

2. **Circos Plot Generation**
After obtaining the `.gbff` file of the new assembly with chromosome description, use `blast_gbff_plot_genes_in_chr.ipynb` to map and visualise the genes in a circos plot.
- Open `blast_gbff_plot_genes_in_chr.ipynb` in Jupyter Notebook
- Follow the notebook instructions to generate circular genome plots
- Output will be saved as PDF files (e.g., `RKQQC_circos.pdf`)

### FindHaplotypesIsoforms
Tools for identifying haplotypes and isoforms in gene sequences.
This uses the consensus sequence of several isolates, following variant calling analysis.

1. **Haplotype Identification**
- Open `gene_haplotype_identify.ipynb` in Jupyter Notebook
- Follow the notebook instructions to analyze gene sequences and identify haplotypes

### PrimerMSAplot
Tools for multiple sequence alignment visualisation.

1. **MSA Visualization**
- Open `fung-viz-msa.ipynb` in Jupyter Notebook
- Follow the notebook instructions to create multiple sequence alignment plots

### TreeBuilding
Scripts for phylogenetic tree construction and sequence analysis.

1. **Tree Building Pipeline**
```bash
sbatch tree_pipeline.sh -i <consensus_directory> -s <minimum percentage of known bases in a sequence for acceptance; suggested: 80> -l <minimum number of accepted samples percentage; suggested: 80>
```

### UpsetPlot
Tools for creating UpSet plots to visualize set intersections.

1. **UpSet Plot Generation**
- Open `upset.ipynb` in Jupyter Notebook
- Load your data in CSV format (see `mock_data.csv` for example)
- Follow the notebook instructions to generate UpSet plots
- Output will be saved as SVG files (e.g., `upsetplot.svg`)

### SIZ1_domains
Analysis tools for SIZ1 protein domains.

1. **Domain Analysis**
- Navigate to the SIZ1_domains directory
- Follow the specific instructions in the directory for domain analysis

### fix-marple-genes
Scripts for fixing gene annotations.

1. **Gene Annotation Fixing**
- Navigate to the fix-marple-genes directory
- Follow the specific instructions in the directory for fixing gene annotations

### search_for_old_genes
Tools for legacy gene identification. For instructions see `scripts/instructions.md`.

1. **Legacy Gene Search**
- Navigate to the search_for_old_genes directory
- Follow the specific instructions in the directory for searching legacy genes 

## Contributing

Feel free to submit issues, fork the repository, and create pull requests for any improvements.

## License

This project is open-source and available for use by the bioinformatics community.

## Note

Please ensure you have all required dependencies installed before running the scripts. Some scripts may require specific input file formats or additional configuration. Check the individual script documentation for specific requirements.
