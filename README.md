# Overview
Identification of alleles was based on three steps，as illustrated in Figure: (1) Identification of paralogs within each haplotype genome;(2)Alleles and non-syntenic homologous gene identification;(3)Other unclassified gene clustering and single gene identification
<img width="883" alt="image" src="https://github.com/user-attachments/assets/3c33a694-e226-4c42-8dff-ff3dcbf8ca18" />

# Requirements
JCVI
Broccoli
diamond
Python version 3.6 or above
numpy
scipy
tqdm


# Usage

1. PolyHomology for allele inditification

        python3 PolyHomology.py -i <input_prefix_file> -p <protein_fasta_file> -c <cds_fasta_file> -g <gff3_file> -ident <identity_percentage> -cov <cover_percentage> -t <thread_count>
         -i / --input: (Required) Path to the prefix file of different subgenomes.
         -p / --protein: (Required) Path to the protein FASTA file.
         -c / --cds: (Required) Path to the CDS FASTA file.
         -g / --gff3: (Required) Path to the GFF3 file.
         -ident / --identity: (Optional, default: 80) Minimum identity percentage to report an alignment.
         -cov / --cover: (Optional, default: 60) Minimum coverage percentage to report an alignment.
         -t / --thread: (Optional, default: 10) Number of threads to use for parallel processing.
   
2. Unclassified gene clustering by Broccoli

        For genes that have not been classified, we use Broccoli (available at https://github.com/rderelle/Broccoli/) to perform further clustering.Broccoli is designed to infer with high precision orthologous groups and pairs of proteins using a mixed phylogeny-network approach.

3. Result

        Syntenic allelic or homoeologous genes ----- allele.blocks.merged 
        Paralogous genes ----- paralogs.txt
        partially non-syntenic allelic or homoeologous genes ----- PolyHomology.result
        use Broccoli to perform further clustering 'PolyUnallelegene.result' file genes to get non-syntenic allelic or homoeologous genes and haplotype-specific genes

# Reference 

        1、Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
        2、Tang, Haibao, Vivek Krishnakumar, Xiaofei Zeng, Zhougeng Xu, Adam Taranto, Johnathan S. Lomas, Yixing Zhang, et al. 2024. “ JCVI: A Versatile Toolkit for Comparative Genomics Analysis.” iMeta 3, e211. https://doi.org/10.1002/imt2.211
        3、Romain Derelle, Hervé Philippe, John K Colbourne, Broccoli: Combining Phylogenetic and Network Analyses for Orthology Assignment, Molecular Biology and Evolution, Volume 37, Issue 11, November 2020,  Pages 3389–3396, https://doi.org/10.1093/molbev/msaa159
   
