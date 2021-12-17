#!/bin/bash

display_usage() { 
	echo -e "\n****   Use this bash script to run MutaPipe   ****\n" 
	echo "IMPORTANT: all three arguments are required and have to be written in double quotes."
	echo -e "\nUsage: $0 [-g GENES, -rsl RELATIVE_SEQUENCE_LENGTH, -cov HSP_COVERAGE] \n-g, --genes 				  Specify genes of interest \n-rsl --relative_sequence_length   	  Set to filter out sequences shorter than a given percentage of the reference sequence (0.1-1.0) \n-cov, --hsp_coverage   			  Set to filter out sequences whose best hsp covers less than a given percentage of the reference sequence (0.1-1.0)\n" 
	echo -e "Usage example: To run MutaPipe for gene NEK1, filtering out sequences covering less than 50% of the reference sequence as well as sequences whose best hsp covers less than 10% of the reference sequence, use: \n$0 \"-g NEK1\" \"-rsl 0.5\" \"-cov 0.1\"\n"
	echo -e "Notes: \n1. To to pass a file containing all gene names use \"-g \$(cat filename)\" as the first argument.\n2. More options are available when running the MutaPipe scripts manually. Use command 'python3 script_name -h' for instructions.\n"
	} 


	# if less than two arguments supplied, display usage 
	if [  $# -le 2 ] 
	then 
		display_usage
		exit 1
	fi 
 
	# check whether user had supplied -h or --help . If yes display usage 
	if [[ ( $# == "--help") ||  $# == "-h" ]] 
	then 
		display_usage
		exit 0
	fi 

	# # this doesn't work for some reason I do not understand...?
	# # check whether user had supplied -g
	# if [[ ( $# == "--genes") ||  $# == "-g" ]]
	# then
	# 	echo "genes are present"
	# 	exit 0
	# fi

echo "Initiating MutaPipe"
echo "==================="

echo "	Input genes: $1"
echo ""
echo "	Sequence length threshold: $2"
echo ""
echo "	HSP coverage threshold: $3"
echo ""

# change to directory where this script and the python scripts are stored
cd "$(dirname "$0")"


# run all MutaPipe python scripts one after another
python3 00_search_pdb.py $1
python3 01_download_files.py
python3 02_parse_cif_files.py
python3 03_parse_fasta_files.py
python3 04_blast_against_reference.py
python3 05_pdb_extract_unsolved_res.py  
python3 06_best_structure_per_mutation.py $2 $3
python3 07_a_ClinVar_Annotations_edirect_per_gene_download_files.py  
python3 07_b_ClinVar_Annotations_edirect_per_gene_parse_files.py 
python3 08_add_clinvar_annotations_to_best_structures.py


#implemented arguments in bash script
# $1 = -g, --genes
# $2 = -rsl, --relative_sequence_length
# $3 = -cov, --hsp_coverage



# Overview arguments for all python scripts (can't be used with this bash script!)
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


