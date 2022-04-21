#!/bin/bash

############################################################
# Set up                                                   #
############################################################

# Set variables
# Get path to the directory where all the python scripts are stored:
MUTAPIPE_DIRECTORY="$(dirname "$0")"

# for some reason it seems we need to store the full path of the MutaPipe directory too
# otherwise we cannot access the uniprot reference sequences later on.
# however, we still need the first variable to read in the default genes (this doesn't work with the full path/second variable)
MUTAPIPE_DIRECTORY_FULL_PATH=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# Default behaviour
GENES="$(cat $MUTAPIPE_DIRECTORY/genes.txt)"
LOG="False"
TARGET_DIRECTORY=$(pwd)
ORGANISM="Homo sapiens"
ALL_PDB_IDS="True"
FORMAT="cif pdb fasta"
POLYPEPTIDES="True"
BLASTp_PATH="blastp"
UNIPROT_REFSEQS="$MUTAPIPE_DIRECTORY_FULL_PATH"/Uniprot_reference_seqs/UP000005640_9606.fasta
RELATIVE_SEQUENCE_LENGTH="0.5"
HSP_COVERAGE="0.1"
N_BEST_STRUCTURES="1"

############################################################
# Help                                                     #
############################################################

Help()
{
   # Display Help
   echo
   echo "****   Use this bash script to run MutaPipe   ****"
   echo
   echo "Usage: $0 [-h|l|t|g|o|a|f|p|b|u|r|c|n]"
   echo
   echo "General options:"
   echo "	-h	HELP				Print this help message and exit"
   echo "	-l	LOG				Write console output to log file in target directory if set to True. Default = $LOG"
   echo "	-t	TARGET_DIRECTORY  		Specify target directory where output will be stored. Default = $TARGET_DIRECTORY"
   echo
   echo "MutaPipe options:"
   echo "	-g	GENES				Specify genes of interest. To to pass a file containing all gene names use -g \"\$(cat filename)\". Default = $GENES"
   echo "	-o	ORGANISM			Set species for which to search pdb structures. Default = $ORGANISM"
   echo "	-a	ALL_PDB_IDS		    	Specify whether to retrieve all (True) or max. 10 PDB IDs (False) per gene. Default = $ALL_PDB_IDS"
   echo "	-f	FORMAT			  	Specify file formats to download. Default = $FORMAT. Options = [cif pdb fasta]"
   echo "	-p	POLYPEPTIDES			Specify whether to extract polypeptide sequence (True) or not (False). Default = $POLYPEPTIDES"
   echo "	-b	BLASTp_PATH			Set path to blastp on your system. Default = $BLASTp_PATH"
   echo "	-u	UNIPROT_REFSEQS			Set path to reference proteome fasta file. Default = $UNIPROT_REFSEQS"
   echo "	-r	RELATIVE_SEQUENCE_LENGTH	Set to filter out sequences shorter than a given % of the reference sequence (0.1-1.0). Default = $RELATIVE_SEQUENCE_LENGTH"
   echo "	-c	HSP_COVERAGE			Set to filter out sequences whose best hsp is shorter than a given % of the reference sequence (0.1-1.0). Default = $HSP_COVERAGE"
   echo "        -n      N_BEST_STRUCTURES               Set number of best structures to be listed in output for each sequence/variant. Default = $N_BEST_STRUCTURES"   
   echo
   echo "Usage examples:"
   echo "(1)"
   echo "To run MutaPipe for gene NEK1, filtering out sequences covering less than 50% of the reference sequence as well as sequences whose best hsp covers less than 10% of the reference sequence, use:"
   echo "$0 -g NEK1 -r 0.5 -c 0.1"	
   echo 
   echo "(2)"
   echo "To run MutaPipe for multiple genes at ones pass a textfile using -g \"\$(cat filename)\" or use double quotes when listing the genes, e.g.:"
   echo "$0 -g \"NEK1 SOD1 FUS\" -r 0.5 -c 0.1"	
   echo
}

############################################################
# Process the input options. Add options as needed.        #
############################################################

# Get the options
while getopts ":hl:t:g:o:a:f:p:b:u:r:c:n:" option; do
   case "${option}" in
      h) # display Help
         Help
         exit;;
      l) # Specify whether to write output to log file
		   LOG=${OPTARG};;
	   t) # Specify target directory 
		   TARGET_DIRECTORY=${OPTARG};;
      g) # Specify input genes
         GENES=${OPTARG};;
	   o) # Set organism of interest 
		   ORGANISM=${OPTARG};;
	   a) # Specify whether to download all or max. 10 PDB IDs per gene
		   ALL_PDB_IDS=${OPTARG};;
	   f) # Specify format of file to be downloaded from the PDB
		   FORMAT=${OPTARG};;
	   p) # Specify whether to extract polypeptide sequences or not
		   POLYPEPTIDES=${OPTARG};;
	   b) # Set PATH to blastp on your system
		   BLASTp_PATH=${OPTARG};;
	   u) # Set PATH to UniProt reference proteome fasta file
		   UNIPROT_REFSEQS=${OPTARG};;
      r) # Set relative sequence lenght
         RELATIVE_SEQUENCE_LENGTH=${OPTARG};;
      c) # Set HSP coverage
         HSP_COVERAGE=${OPTARG};;
      n) # Set number of best structures per gene/sequence to be returned in the final output table
         N_BEST_STRUCTURES=${OPTARG};;
     \?) # Invalid option
         echo "Error: Invalid option."
         echo "Use -h to display help message."
         exit;;
   esac
done

############################################################
############################################################
# MutaPipe 		                                           #
############################################################
############################################################

echo
echo "Initiating MutaPipe"
echo "==================="
echo
echo "Input parameters:	GENES 				$GENES"
echo "			LOG 				$LOG"
echo "			TARGET_DIRECTORY 		$TARGET_DIRECTORY"
echo "			ORGANISM 			$ORGANISM"
echo "			ALL_PDB_IDS			$ALL_PDB_IDS"
echo "			FORMAT				$FORMAT"
echo "			POLYPEPTIDES			$POLYPEPTIDES"
echo "			BLASTp_PATH			$BLASTp_PATH"
echo "			UNIPROT_REFSEQS 		$UNIPROT_REFSEQS"
echo "			RELATIVE_SEQUENCE_LENGTH	$RELATIVE_SEQUENCE_LENGTH"
echo "			HSP_COVERAGE: 			$HSP_COVERAGE"
echo "                        N_BEST_STRUCTURES:              $N_BEST_STRUCTURES"


# change to directory where this script and the python scripts are stored
cd "$MUTAPIPE_DIRECTORY"

# run all MutaPipe python scripts one after another
python3 00_search_pdb.py -g $GENES -o "$ORGANISM" -a $ALL_PDB_IDS -t "$TARGET_DIRECTORY" -l $LOG 
python3 01_download_files.py -f $FORMAT -t "$TARGET_DIRECTORY" -l $LOG
python3 02_parse_cif_files.py -pp $POLYPEPTIDES -t "$TARGET_DIRECTORY" -l $LOG
python3 03_parse_pdb_files_extract_unsolved_residues.py -t "$TARGET_DIRECTORY" -l $LOG
python3 04_parse_fasta_files.py -t "$TARGET_DIRECTORY" -l $LOG
python3 05_blast_against_reference.py -bp "$BLASTp_PATH" -refseq "$UNIPROT_REFSEQS" -t "$TARGET_DIRECTORY" -l $LOG
python3 06_a_download_ClinVar_data.py -t "$TARGET_DIRECTORY" -l $LOG
python3 06_b_parse_ClinVar_data.py -t "$TARGET_DIRECTORY" -l $LOG
python3 07_combine_data_to_get_best_n_structures_per_sequence.py -rsl $RELATIVE_SEQUENCE_LENGTH -cov $HSP_COVERAGE -t "$TARGET_DIRECTORY" -l $LOG -n_best $N_BEST_STRUCTURES


# change (back) to target directory
cd "$TARGET_DIRECTORY"


############################################################
# Overview arguments for all python scripts			       #
############################################################

# options for all scripts
#   -h, --help                      show help message and exit
#   -l, --log 			               Write console output to log file in current directory if set to True, default = False
#   -t, --target TARGET             Specify target directory, default = /Users/debs/OneDrive - King's College London/Pipeline/Pipeline_Git/MutaPipe

# additional options for script 00_search_pdb.py
#   -g, --genes 		               Specify genes for which to search pdb structures, default = ['OPTN', 'ERBB4', 'DCTN1']; to pass a file containing all genes use -g $(cat filename)
#   -o, --organism		            Specify species for which to search pdb structures, default = Homo sapiens
#   -a, --all		                  Retrieve all (True) vs max. 10 pdb IDs per gene (False), default = True

# additional options for script 01_download_files.py
#   -f, --format			            Specify file format to be downloaded. For mmCif files (.cif) use 'cif' ; for pdb files (.pdb) use 'pdb' ; for fasta files (.fasta) use 'fasta' ; default = cif pdb fasta

# additional options for script 02_parse_cif_files.py
#   -pp, --polypeptides             Specify whether to extract polypeptide sequence (True) or not (False), default = True

# new in 05_blast_against_reference.py
# -bp, --blastp_path 	            Specify the path to blastp on your system ; default = /usr/local/ncbi/blast/bin/blastp       
# -refseq, --reference_sequences    Specify path to uniprot reference fasta, default = /Users/debs/OneDrive - King's College London/Pipeline/Pipeline_Git/MutaPipe/Uniprot_reference_seqs/UP000005640_9606.fasta   				                   

# new in 07_combine_data_to_get_best_n_structures_per_sequence.py
# -rsl, --relative_sequence_length  filter out sequences shorter than a given percentage of the reference sequence, default = 0.5
# -cov, --hsp_coverage              filter out sequences whose best hsp covers less than a given percentage of the reference sequence, default = 0.1
# -n_best, --n_best_structures      number of best structures to be listed in output for each sequence/variant, default = 1



