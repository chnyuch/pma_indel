# pma_indel
The pipeline for indel identification 

1. Download the gff files and the genome fasta files of the focal species
2. Use AGAT to acquire the longest isoform of the cds and the protein fasta
3. orthoFinder_cmd.1.py: Use Orthofinder to acquire the orthologous groups between the focal species
4. orthoFinder_prep.2.py: Organize the output from orthoFinder
5. indel_aln_recon.2.py: Reconstruct ancestral state of the sequences using indel-aware methods (indelMap)
6. indel_iden.3.py: Identify the indels and filter out the great tit specific indel
7. indel_selPrep.5.py: Acquire the new alignments based on the sites to the indels and use PAML to calculate dnds.
8. indel_parsePaml.1.py: parsing the paml output
