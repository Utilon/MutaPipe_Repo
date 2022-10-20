# This script takes a csv files as input, namely:
#      - 01_search_overview_folders.csv
# and will
#      - search for pdb files in all folders listed in 01_search_overview_folders.csv
#      - extract info on missing residues / residues which have not been solved in the crystal structure for each pdb structure
#      - output the following files:
#                - a csv file called 03_unsolved_residues_per_structure.csv listing all unsolved residues in all structures for all genes
#                  (one row for each structure)
#                - a csv files called 03_unsolved_residues_per_chain.csv listing all unsolved residues in all chains of all structures for all genes
#                  (one row for each chain)
#  ----------------------------------------------------------------------------------------------------------------------------------
   
# Set up
import pandas as pd
import os
from os import listdir
from os.path import isfile, join, exists
import ast
import sys
import argparse
from datetime import datetime

from Bio.PDB import *

# get this script's name:
script_name = os.path.basename(__file__)

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
target_directory = os.getcwd()    # set target directory (where Results folder is located)
delete_files=False                  # specify whether to delete pdb files after parsing them or not
# download_files_to_separate_directory = True # specify if pdb mmcif and fasta files should be stored in separate directory
web_run = True # specify if pdb mmcif and fasta files should be stored in separate directory
mutafy_directory = f'{target_directory}/mutafy' # set path to folder where structures will be/are stored
                                                                                        
                                            
# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****    This script takes a csv files as input, namely:
01_search_overview_folders.csv and will:
1. search for pdb files in all folders listed in 01_search_overview_folders.csv
2. extract info on missing residues / residues which have not been solved in the crystal structure for each pdb structure
4. output the following files:
(1) a csv file called 03_unsolved_residues_per_structure.csv listing all unsolved residues in all structures for all genes (one row for each structure)
(2) a csv files called 03_unsolved_residues_per_chain.csv listing all unsolved residues in all chains of all structures for all genes (one row for each chain)    ***""")

ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')
ap.add_argument("-del", "--delete_files", type=str2bool, required = False, help=f'Specify whether to delete pdb files after parsing (True) or not (False), default = {str(delete_files)}')
# ap.add_argument("-s", "--sep_dir", type=str2bool, required = False, help=f'specify if pdb mmcif and fasta files should be stored in separate directory (True) or not (False), default = {str(download_files_to_separate_directory)}')
ap.add_argument("-w", "--web_run", type=str2bool, required = False, help=f'Indicate whether MutaPipe is run via a webserver (True) or not (False), default = {str(mutafy_directory)}')
ap.add_argument("-m", "--mutafy", required = False, help=f'set path to mutafy directory where information from previous runs is stored, default = {mutafy_directory}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
delete_files  = delete_files if args["delete_files"]   == None else args["delete_files"]
web_run = web_run if args["web_run"] == None else args["web_run"]
mutafy_directory = mutafy_directory if args["mutafy"] == None else args["mutafy"]

# ----------------------------------------------------------------------------------------------------------------------------------
# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# ----------------------------------------------------------------------------------------------------------------------------------
#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_03.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_03.txt', 'a')

# print nice title
print('===============================================================================')
print('*****    Extracting Unsolved Residues from PDB Structures for Input Genes    *****')
print('===============================================================================\n')

# print script name to console/log file
print(f'script name: {script_name}')


# store current date and time in an object and print to console / write to log file
start_time = datetime.now()
print(f'start: {start_time}\n')

# ----------------------------------------------------------------------------------------------------------------------------------

# if there is no data in the PDB for any of the input genes, we don't have to run this script (or any of the
# following scripts in the pipeline, apart from the AlphaFold one)
# so we read in the file 00_search_overview_availability.csv from the results_dir to check if we have to run the script
data_availability = pd.read_csv(f'{results_dir}/00_search_overview_availability.csv')
# the data_availability df has two columns, one with the gene_name and one with a boolean value indicating if PDB data is available for this gene
# if all values in the data_available column are False, we can skip this script
if True not in data_availability.data_available.unique():
    print('No PDB data available for any of the input genes')
    print ('Exiting Python...')
    sys.exit('No PDB data available for any of the input genes')


# read in csv files
folders = pd.read_csv(f'{results_dir}/01_search_overview_folders.csv')

# if this is a web run, we also read in the csv file from script 01 which lists all the new pdb ID's to be downloaded / parsed!
if web_run:
    df_new_structures_to_download = pd.read_csv(f'{mutafy_directory}/01_new structures_to_be parsed_mutafy.csv')
    
# create an empty df called unsolved to populate with information on unsolved residues extracted from the pdb files in each folder
unsolved = pd.DataFrame(columns=['gene', 'structure_id', 'unsolved_residues_in_structure'])

# we also create a df to store similar information, but here we want one row for each chain, instead of one row for each structure
unsolved_per_chain = pd.DataFrame(columns=['gene', 'structure_id', 'chain', 'unsolved_residues_in_chain'])

# define variable to keep track of genes being looked at (for console output)
previous_gene = 'no_gene'
counter = 0

# define variable to count number of pdb files parsed overall
pdb_total = 0

# we loop over the folders df to check each of the listed folders for pdb files
for index, row in folders.iterrows():
    counter += 1
    gene = row.gene_name
    structure_folder = row.full_path
    # change to relevant folder
    os.chdir(structure_folder)
    
    # if this is a webrun, we first check if there are any new structures to be parsed in this folder/for this gene:
    if web_run:
        # get the list of pdb ids to be parsed        
        new_structures_to_download = ast.literal_eval(df_new_structures_to_download[df_new_structures_to_download.gene == gene].new_pdb_ids.values[0])
        # currently this list contains pdb ids, but in order for the rest of the loop to work, we need a variable called
        # pdb_files which contains pdb filenames to be parsed (pdb.pdb)
        # so we do the following:
        pdb_files = [f'{pdb_id}.pdb' for pdb_id in new_structures_to_download]
        # if there are no new structures to be parsed, we can continue to the next gene/folder
        if len(pdb_files) == 0:
            print(f'\nNo new pdb files to be parsed for {gene} (gene {counter} of {len(folders)})\n')
            continue
        
    elif web_run == False:      
#         # change to the folder containing the files to be parsed
#         os.chdir(structure_folder)              
        # create list with filenames of all pdb/mmCIF files in this folder
        files = [f for f in listdir(row.full_path) if isfile(join(row.full_path, f))]
        pdb_files = [f for f in files if '.pdb' in f]
        # sort list
        pdb_files.sort()
            
    # update pdb_total
    pdb_total += len(pdb_files)
    
    # print information to console
    if gene != previous_gene:
        print(f'>>> Looping over {len(pdb_files)} pdb files for gene {gene} (gene {counter} of {len(folders)})...')
    
    # now that we have a list of all pdb files in the current folder, we can extract information from them and append it to our unsolved df
    
    # Biopython has a method to parse pdb headers
    # The missing residues are also stored in the pdb header information
    print('    >>> parsing pdb headers and extracting unsolved residues...\n')
    
    # as we want to do this for all our pdb files, we initiate a for loop:
    for pdb in pdb_files:
        # added try and except statement as sometimes the pdb files for newly available structures in a webrun
        # are not available, e.g. no pdb file because structure is too large.
        # if it's not a webrun, then this is not a problem, because we parse all pdb files in a given folder, but in case of a webrun,
        # we specifiy the new ones, so that's why.
        try:
            # we parse the header like so:
            header = parse_pdb_header(pdb)
        except FileNotFoundError:
            print(f'No file {pdb} exists. This is likely due to the structure being too large for pdb file format.')
            # substract -1 from the pdb_total:
            pdb_total -= 1
            continue
        # now that we have parese the header information from the pdb file, we can delete the pdb file (don't use it anymore)
        if delete_files == True:
            os.remove(pdb)
        # and we use the header to get the missing residues
        missing_res = header['missing_residues']
        # missing_res is a list of dictionaries. For each unsolved/missing residue, there is one dictionary with the following keys:
        # 'model', 'res_name', 'chain', 'sseq', 'insertion'
        
        # create empty dictionary to store unsolved residues and positions per chain
        missing_res_dict = {}
        
        # loop over missing_res dictionary to fill missing_res_dict dictionary
        # for each residue we get the name, chain and position
        for residue in missing_res:
            name = residue['res_name']
            chain = residue['chain']
            pos = residue['ssseq']
            if chain in missing_res_dict.keys():
                # if the chain is a key in missng_res_dict, we append this residue to the corresponding value (list)
                missing_res_dict[chain].append(f'{name}{pos}')
            elif chain not in missing_res_dict.keys():
                # if the chain is not a key, we create the key first and then add the value (in list format)
                missing_res_dict[chain] = [f'{name}{pos}']
            
        # now that we have a dictionary with all the missing residues per chain for the current structure,
        # we can append a row to the df 'unsolved'
        # columns of unsolved are ['gene', 'structure_id', 'unsolved_residues_in_structure']
        unsolved.loc[len(unsolved)] = [gene, pdb[:4], missing_res_dict]
        
        # we also make a df containing one row per chain (insted of one row per structure), which we can later combine more easily with the info from the blastp output (which are also per chain)
        for key, value in missing_res_dict.items():
                unsolved_per_chain.loc[len(unsolved_per_chain)] = [gene, pdb[:4], key, value]

    if not web_run:
        # now that we have looped over all the pdb files in this folder and
        # extracted the unsolved residues and saved them in the unsolved_per_chain and the unsolved df,
        # we can extract slices from the current versions of the unsolved_per_chain and unsolved df for the current gene
        unsolved_this_gene = unsolved[unsolved.gene == gene]
        unsolved_per_chain_this_gene = unsolved_per_chain[unsolved_per_chain.gene == gene]
        # we can save this to the current folder (gene folder)
        unsolved_this_gene.to_csv(f'{gene}_03_unsolved_residues_per_structure.csv', index=False)
        unsolved_per_chain_this_gene.to_csv(f'{gene}_03_unsolved_per_chain.csv', index=False)

# change to results directory
os.chdir(results_dir)

# now, if this is a webrun, we want to read in the data from previous mutafy runs and
# get a slice for all the genes of the current seach (otherwise only new structures will be listed in the output file)
# also, we want to update already existing mutafy data and add new structures to it!
if web_run:
    if exists(f'{mutafy_directory}/03_unsolved_residues_per_structure_mutafy.csv'):
        mutafy_unsolved_res = pd.read_csv(f'{mutafy_directory}/03_unsolved_residues_per_structure_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        # in order to do that, we have to convert the values in the unsolved df in column unsolved_residues_in_structure
        # to a string (currently this is in dict format)
        unsolved["unsolved_residues_in_structure"]=unsolved["unsolved_residues_in_structure"].values.astype('str')
        # now we can concatenate the two dfs
        updated_mutafy_unsolved_res = pd.concat([mutafy_unsolved_res, unsolved], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name and then structure id
        updated_mutafy_unsolved_res.sort_values(by=['gene', 'structure_id'], inplace=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_unsolved_res.to_csv(f'{mutafy_directory}/03_unsolved_residues_per_structure_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the resolutions for all the structures for all genes of the current webrun
        unsolved = updated_mutafy_unsolved_res[updated_mutafy_unsolved_res.gene.isin(list(folders.gene_name.values))]
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        unsolved.to_csv(f'{mutafy_directory}/03_unsolved_residues_per_structure_mutafy.csv', index = False)

    if exists(f'{mutafy_directory}/03_unsolved_residues_per_chain_mutafy.csv'):
        mutafy_unsolved_res_per_chain = pd.read_csv(f'{mutafy_directory}/03_unsolved_residues_per_chain_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        # in order to do that, we convert all the values in the unsolved_per_chain df to strings
        unsolved_per_chain = unsolved_per_chain.astype('str')
        updated_mutafy_unsolved_res_per_chain = pd.concat([mutafy_unsolved_res_per_chain, unsolved_per_chain], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name and then structure id
        updated_mutafy_unsolved_res_per_chain.sort_values(by=['gene', 'structure_id'], inplace=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_unsolved_res_per_chain.to_csv(f'{mutafy_directory}/03_unsolved_residues_per_chain_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the resolutions for all the structures for all genes of the current webrun
        unsolved_per_chain = updated_mutafy_unsolved_res_per_chain[updated_mutafy_unsolved_res_per_chain.gene.isin(list(folders.gene_name.values))]
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        unsolved_per_chain.to_csv(f'{mutafy_directory}/03_unsolved_residues_per_chain_mutafy.csv', index = False)
        
# write output to csv files
print('>>> writing csv file containing unsolved residues per structure for all genes...')
print('>>> writing csv file containing unsolved residues per chain for all structures for all genes...\n')
unsolved.to_csv('03_unsolved_residues_per_structure.csv', index=False)
unsolved_per_chain.to_csv('03_unsolved_residues_per_chain.csv', index=False)

# change back to target_directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Parsed a total of {pdb_total} pdb files stored across {len(folders)} folders.')

print('\nThe following files have been created and stored in the Results folder:')
print('   o      03_unsolved_residues_per_structure.csv     (lists unsolved residues per structure for all genes (one row = one structure))')
print('   o      03_unsolved_residues_per_chain.csv           (lists unsolved residues per chain for all genes (one row = one chain))\n')

# print script name to console/log file
print(f'end of script {script_name}')

# store current date and time in an object and print to console / write to log file
end_time = datetime.now()
print(f'start: {start_time}')
print(f'end: {end_time}\n\n')
print('........................................................................................................................................................\n\n\n')

# close search log
if create_search_log == True:
    sys.stdout.close()