# Set up

import sys
import argparse
from datetime import datetime

# ----------------------------------------------------------------------------------------------------------------------------------

# use argparse to make it so we can pass arguments to script via terminal

# define a function to convert different inputs to booleans
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# set default values for arguments we want to implement
# we have to do this here if we want to print the default values in the help message

create_search_log = False     # will create a file called search_log.txt with console output if set to True,
                                            # prints to console if set to False.
download_format =  'cif pdb fasta'      # specify which file formats to download from the PDB for the input data
                                            # use the following abbreviations 
                                                # for mmCif files (.cif) use 'cif'
                                                # for pdb files (.pdb) use 'pdb'
                                                # for fasta files (.fasta) use 'fasta'

target_directory = os.getcwd()    # set target directory (where Results folder is located)
                                            
                                            
# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****    This script takes a csv file containing gene names and corresponding PDB IDs as input and will:
1. create a folder for each gene in the results folder
2. download specified formats (cif, pdb, fasta) into respective folders (atm automatically downloads mmCIF, pdb and FASTA files)
3. outputs a csv file called 01_search_overview_folders listing all the the newly created folders and their contents 
4. outputs a csv file called 01_search_overview_n_structures.csv listing the number of structures retrieved per gene    ***""")

ap.add_argument('-f','--format', nargs='+', required=False, help=f"Specify file format to be downloaded. For mmCif files (.cif) use 'cif' ; for pdb files (.pdb) use 'pdb' ; for fasta files (.fasta) use 'fasta' ; default = {download_format}")
ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in current directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
download_format = download_format if args["format"] == None else args["format"]
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]

# ----------------------------------------------------------------------------------------------------------------------------------

if 'cif' in download_format:
    print('downloading mmCif...')
if 'fasta' in download_format:
    print('downloading fasta...')
if 'pdb' in download_format:
    print('downloading pdb...')