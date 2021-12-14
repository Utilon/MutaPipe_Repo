#!/bin/bash


python3 00_search_pdb.py $1
python 01_download_files.py
python 02_parse_cif_files.py
python 03_parse_fasta_files.py
python 04_blast_against_reference.py
python 05_pdb_extract_unsolved_res.py  
python 06_best_structure_per_mutation.py 
python 07_a_ClinVar_Annotations_edirect_per_gene_download_files.py  
python 07_b_ClinVar_Annotations_edirect_per_gene_parse_files.py 
python 07_add_clinvar_annotations_to_best_structures.py




#implemented arguments in bash script
# $1 = -g, --genes
# $2 = -rsl, --relative_sequence_length
# $3 = -cov, --hsp_coverage



# Overview arguments for all python scripts
#==========================================
# options for all scripts
#   -h, --help            show help message and exit
#   -l, --log 			  Write console output to log file in current directory if set to True, default = False
#   -t, --target TARGET   Specify target directory, default = /Users/debs/OneDrive - King's College London/Pipeline/Pipeline_Git/MutaPipe


# options for script 00_search_pdb.py
#   -g, --genes 		  Specify genes for which to search pdb structures, default = ['OPTN', 'ERBB4', 'DCTN1']; to pass a file containing all genes use -g $(cat filename)
#   -o, --organism		  Specify species for which to search pdb structures, default = Homo sapiens
#   -a, --all		      Retrieve all (True) vs max. 10 pdb IDs per gene (False), default = True


# new in 01
#  -f, --format			  Specify file format to be downloaded. For mmCif files (.cif) use 'cif' ; for pdb files (.pdb) use 'pdb' ; for fasta files (.fasta) use 'fasta' ; default = cif pdb fasta


# new in 02
# -pp, --polypeptides 	  Specify whether to extract polypeptide sequence (True) or not (False), default = True


# new in 04
# -bp, --blastp_path 	  Specify the path to blastp on your system ; default = /usr/local/ncbi/blast/bin/blastp       
# -refseq, --reference_sequences
#   				      Specify path to uniprot reference fasta, default = /Users/debs/OneDrive - King's College London/Pipeline/Pipeline_Git/MutaPipe/Uniprot_reference_seqs/UP000005640_9606.fasta    

# new in 06
# -rsl, --relative_sequence_length
#                         filter out sequences shorter than a given percentage of the reference sequence, default = 0.5
# -cov, --hsp_coverage    filter out sequences whose best hsp covers less than a given percentage of the reference sequence, default = 0.1


