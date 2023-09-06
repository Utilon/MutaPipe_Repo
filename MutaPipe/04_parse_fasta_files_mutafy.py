# This script takes a csv file (01_search_overview_folders.csv) containing information on the folders
# where the fasta files and fasta_ex files are stored as input and will:
#      - extract information from each fasta file, incl.:
#                - chain name
#                - description
#                - species 
#                - sequence
#      - extract information from each fasta_ex file, incl.:
#                - chain name
#                - description
#                - uniprot id 
#                - sequence
#      - merge the data and combine information from both fasta and fasta_ex files
#      - output the following files:
#                - a csv file called 04_fasta_info.csv containing information extracted from all fasta files for all genes
#                - a csv file called 04_fasta_ex_info.csv containing information extracted from all fasta_ex files for all genes
#                - a csv file called 04_fasta_combined_info.csv containing combined information extracted from all fasta and fasta_ex files for all genes    
#  ----------------------------------------------------------------------------------------------------------------------------------
   
# Set up
import pandas as pd
import os
from os import listdir
from os.path import isfile, join, exists
import ast

from Bio import SeqIO
from Bio.PDB import *
import sys
import argparse
from datetime import datetime

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
delete_files=False                  # specify whether to delete fasta files after parsing them or not
# download_files_to_separate_directory = True # specify if pdb mmcif and fasta files should be stored in separate directory
web_run = True # specify if pdb mmcif and fasta files should be stored in separate directory
mutafy_directory = f'{target_directory}/mutafy' # set path to folder where structures will be/are stored

# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****    This script takes a csv file (01_search_overview_folders.csv) containing information on the folders where all fasta files are stored as input and will: 
1. extract information from each fasta file 
2. merge information from multiple fasta files whenever available (fasta and fasta_ex files [extracted from mmCif]) 
3. output the following files: 
(1) a csv file called GENENAME_04_fasta_info.csv per gene/folder containing information extracted from all fasta files for this gene 
(2) a csv file called GENENAME_04_fasta_ex_info.csv per gene/folder containing information extracted from all fasta_ex files for this gene 
(3) a csv file called GENNAME_04_fasta_combined_info.csv per gene/folder containing combined information extracted from all fasta and fasta_ex files for this genes 
(4) a csv file called 04_fasta_info.csv containing information extracted from all fasta files for all genes 
(5) a csv file called 04_fasta_ex_info.csv containing information extracted from all fasta_ex files for all genes 
(6) a csv file called 04_fasta_combined_info.csv containing combined information extracted from all fasta and fasta_ex files for all genes    ***""")

ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')
ap.add_argument("-del", "--delete_files", type=str2bool, required = False, help=f'Specify whether to delete pdb files after parsing (True) or not (False), default = {str(delete_files)}')
ap.add_argument("-w", "--web_run", type=str2bool, required = False, help=f'Indicate whether MutaPipe is run via a webserver (True) or not (False), default = {str(mutafy_directory)}')
ap.add_argument("-m", "--mutafy", required = False, help=f'set path to mutafy directory where information from previous runs is stored, default = {mutafy_directory}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
delete_files = delete_files if args["delete_files"] == None else args["delete_files"]
web_run = web_run if args["web_run"] == None else args["web_run"]
mutafy_directory = mutafy_directory if args["mutafy"] == None else args["mutafy"]
# ----------------------------------------------------------------------------------------------------------------------------------
# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# ----------------------------------------------------------------------------------------------------------------------------------
#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_04.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_04.txt', 'a')

# print nice title
print('===============================================================================')
print('*****    Parsing Fasta Files from the Protein Data Bank for Input Genes    *****')
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

# Read in data from csv file
folder_info = pd.read_csv(f'{results_dir}/01_search_overview_folders.csv')#, usecols=['folder_name', 'full_path'])

# if this is a web run, we also read in the csv file from script 01 which lists all the new pdb ID's to be downloaded / parsed!
if web_run:
    df_new_structures_to_download = pd.read_csv(f'{mutafy_directory}/01_new structures_to_be parsed_mutafy.csv')

# define variable n_folders (number of folders to be checked for files) = number of rows in dataframe
n_folders = len(folder_info)

# create variable to keep track of how many fasta and fasta_ex files get parsed overall
fasta_total = 0
fasta_ex_total = 0

# ultimately we want to create dfs with information from all fasta files and all fasta_ex files
# need to extract the following:
#                                  'chain_id', # ??? chain counter
#                                  'chain_name', # name of chain (extracted from fasta file)
#                                  'species', # name of species associated with this chain/sequence (extracted from fasta file)                                 
#                                  'description', # description (extracted from fasta files)
#                                  'description_code' # description code (extracted from fasta_ex file)
#                                  'uniprot_id', # uniprot id of this sequence (extracted from fasta_ex file)

# make two empty lists to populate with information from all fasta files and all fasta_ex files
fasta_list = []
fasta_ex_list = []
# and one list for the combined info from fasta and fasta_ex files
combined_list = []

# we loop over the df containing the folder names and the full paths to each folder
# we use this information to parse the files in each folder and extract information
folder_counter = 0

for index, row in folder_info.iterrows():
    folder_counter += 1
    gene = row.gene_name 
    structure_folder = row.full_path
    # change to the folder containing the files to be parsed
    os.chdir(structure_folder)   
    
    # create list with filenames of all fasta and fasta_ex files in this folder
    files = [f for f in listdir(structure_folder) if isfile(join(structure_folder, f))]
    fasta_files = [file for file in files if ('.fasta' in file) & ('_ex.' not in file)]
    fasta_ex_files = [file for file in files if ('.fasta' in file) & ('_ex.' in file)]
    
    print(f"""\nChecking folder content (folder {folder_counter} of {n_folders}) for gene {gene}:""")
    print(f"            number of fasta files ('.fasta'):                       {len(fasta_files)}")
    print(f"            number of fasta_ex files ('_ex.fasta'):            {len(fasta_ex_files)}")    
    
    # if this is a webrun, we check if there are any *new* fasta files to be parsed in this folder/for this gene:
    if web_run:
        # get the list of pdb ids to be parsed        
        new_structures_to_download = ast.literal_eval(df_new_structures_to_download[df_new_structures_to_download.gene == gene].new_pdb_ids.values[0])
        # currently this list contains pdb ids, but in order for the rest of the loop to work, we need a variable called
        # pdb_files which contains fasta filenames to be parsed (pdb_id.fasta)
        # so we do the following:
        fasta_files = [f'{pdb_id}.fasta' for pdb_id in new_structures_to_download]
        fasta_ex_files = [f'{pdb_id}_ex.fasta' for pdb_id in new_structures_to_download]
    
    # if there are no new fasta files or fasta_ex files to be parsed, we can print a statement, and continue to the next folder
    if len(fasta_files)==0 & len(fasta_ex_files)==0:
        print(f'\nNo new fasta or fasta_ex files have to be parsed for {gene} (gene {folder_counter} of {len(folder_info)})\n')
        continue
        
    # sort lists
    fasta_files.sort()
    fasta_ex_files.sort()
    
    # add number of fasta and fasta_ex files to be parsed to total number of fasta and fasta_ex files, respectively:
    fasta_total += len(fasta_files)
    fasta_ex_total += len(fasta_ex_files)
    
# THIS SECTION IS INCLUDED ABOVE NOW, DELETE HERE IF EVERYTHING WORKS AS INTENDED!
#     # if there are no new fasta files to be parsed, we can continue to the next gene/folder
#     if len(fasta_files) == 0:
#         print(f'\nNo new fasta or fasta_ex files have to be parsed for {gene} (gene {folder_counter} of {len(folder_info)})\n')
#         continue

# THIS SECTION IS INCLUDED ABOVE NOW, DELETE HERE IF EVERYTHING WORKS AS INTENDED!
#     print(f"""\nChecking folder content (folder {folder_counter} of {n_folders}) for gene {gene}:""")
#     print(f"            number of fasta files ('.fasta'):                       {len(fasta_files)}")
#     print(f"            number of fasta_ex files ('_ex.fasta'):            {len(fasta_ex_files)}")    
    
    # LOOP OVER FASTA FILES
    # -------------------------------------
    # now we can loop over all fasta files in this folder to extract information and append it to the fasta_df
    print(f'\n{len(fasta_files)} fasta file(s) have to be parsed for {gene} (gene {folder_counter} of {len(folder_info)})')
    for fasta in fasta_files:
        print(f'    >>> parsing {fasta}')
        # added try and except statement as sometimes the fasta files for newly available structures in a webrun
        # may not be available (I think).
        # if it's not a webrun, then this is not a problem, because we only parse fasta files in a given folder, but in case of a webrun,
        # we specifiy the new ones, so that's why.
        try:
            # we read in all sequences/chains from a fasta file and store them in fasta_records 
            fasta_records = list(SeqIO.parse(fasta, 'fasta'))
        except FileNotFoundError:
            print(f'No fasta file detected for {fasta[:4]} for gene {gene}!\n')
            continue
        # now that we've read in the sequence from the fasta file, we can delete it (don't use it anymore)
        if delete_files == True:
            os.remove(fasta)
        # now we loop over all sequence records to extract the information associated with every sequence
        record_id = 0
        for record in fasta_records:
            record_id += 1
            # first we extract information from the record description: chain_name, description, species
            # record.description comes in the following format:
            # '6KJ2_1|Chain A|RNA-binding protein FUS|Homo sapiens (9606)'
            chain_name = record.description.split('|')[1]
            description = record.description.split('|')[2]
            species = record.description.split('|')[3]
            # adjust format of description and species if it's in all caps
            if description.isupper() or description.islower():
                description = description.capitalize()
            if species.isupper() or species.islower():
                species = species.capitalize()            
            # we also want to extract the sequence for this record/chain
            sequence = record.seq
            
            # now we add all this data to the fasta_list which will be converted to a df later
            fasta_list.append({'gene_name': gene,
                              'structure_id': fasta.replace('.fasta', ''),
                              'record_id': record_id,
                              'chain_name': chain_name,
                              'species': species,
                              'description': description,
                              'sequence': sequence})
         
    # LOOP OVER FASTA_EX FILES
    # -------------------------------------
    # we do the same for the fasta_ex files:
    # we  loop over all fasta_ex files in this folder to extract information and append it to the fasta_ex_df
    print(f'{len(fasta_ex_files)} fasta_ex file(s) have to be parsed for {gene} (gene {folder_counter} of {len(folder_info)})')
    for fasta_ex in fasta_ex_files:
        print(f'    >>> parsing {fasta_ex}')
        # added try and except statement as sometimes the fasta files for newly available structures in a webrun
        # may not be available.
        try:        
            # we read in all sequences/chains from a fasta_ex file and store them in fasta_ex_records 
            fasta_ex_records = list(SeqIO.parse(fasta_ex, 'fasta'))
        except FileNotFoundError:
            print(f'No fasta_ex file detected for {fasta[:4]} for gene {gene}')
            continue
        # now that we've read in the sequence from the fasta_ex file, we can delete it (don't use it anymore)
        if delete_files == True:
            os.remove(fasta_ex)        
        # now we loop over all sequence records to extract the information associated with every sequence
        record_id_ex = 0
        for record in fasta_ex_records:
            record_id_ex += 1            
            # first we want to extract the sequence for this record/chain
            sequence_ex = record.seq
            # then we extract information from the record description: chain_name, description, uniprot_id
            # record.description comes in the following format:
            # '6KJ2:A UNP:P35637 FUS_HUMAN'
            # sometimes the record.description in fasta_ex files is unknown and comes in the format
            # A <unknown description>
            # if this is the case, we set the following variables to unknown
            if 'unknown' in record.description:
                chain_name_ex = 'unknown' 
                uniprot_id_ex = 'unknown'
                description_ex = 'unknown'
            # also, sometimes the record.description in fasta_ex files is weird and has the following format:
            # 6G99:A PDB:6G99 6G99
            # if this is the case, record.description.split(' ')[1].split(':')[1].lower() should be equal to the current structure id:
            elif record.description.split(' ')[1].split(':')[1].lower() == fasta_ex.replace('_ex.fasta', '').lower():
                # in this case, we also want to set the same variables to unknown:
                chain_name_ex = 'unknown'
                uniprot_id_ex = 'unknown'
                description_ex = 'unknown'                
            # if the description is not unknow or in this weird format, we do the following instead to extract information:
            # in case and record.description comes in another weird format which can't be parsed like this, we add a try except statement:
            else:
                try:
                    chain_name_ex = record.description.split(' ')[0].split(':')[1]  # this is not the actual chain name, but it just starts with 'A' for the first and then goes on alphabetically
                    uniprot_id_ex = record.description.split(' ')[1].split(':')[1]
                    description_ex = record.description.split(' ')[2]
                # if this doesn't work because the description is in some other format we don't capture, we just set the variables to 'unknown'
                except:
                    chain_name_ex = 'unknown' 
                    uniprot_id_ex = 'unknown'
                    description_ex = 'unknown'                
            
            # now we append all the fasta_ex data to our fasta_ex_list
            fasta_ex_list.append({'gene_name': gene,
                                  'structure_id': fasta_ex.replace('_ex.fasta', ''),
                                  'record_id': record_id_ex,
                                  'chain_name': chain_name_ex,
                                  'uniprot_id': uniprot_id_ex,
                                  'description': description_ex,
                                  'sequence': sequence_ex})

# COMBINE FASTA AND FASTA_EX INFORMSTION
# now that we have all data for all genes (fasta and fasta_ex), we convert them to dfs and combine them
fasta_df = pd.DataFrame(fasta_list).astype('str') # convert all values to string for concatenate to work later
fasta_ex_df = pd.DataFrame(fasta_ex_list).astype('str')
# in fasta_ex_df, each chain is listed in a separate row even if it has the same sequence as other chains in the same structure
# in this case, all the information except for the record_id, and the chain_name are the same for all these rows.
# Problem: sometimes the rows are the same as described above (same gene, sequence, structure_id)
# but some might have an unknown  description/uniprot_id while others may not
# In order to keep the rows with valid description/uniprot_id we first have to sort the df so that the row with valid a
# description/uniprot_id always comes first (before the one with unknown description/uniprot_id)
# So we make an extra col with boolean values if uniprot_id and description are 'unknown'
# they are always either both unknown or both not unknown

# We add an extra col in fasta_ex_df to indicate if each row contains unknown stuff
fasta_ex_df['contains_unknowns'] = fasta_ex_df.description.apply(lambda x: True if 'unknown' in x else False)

# now we can sort the df and also use the new column 'contains_unknowns'
# the resulting df will list duplicated rows *without* unknowns first ('False' comes before 'True' in 'contains_unknowns' column)
fasta_ex_df.sort_values(['gene_name', 'structure_id', 'sequence', 'contains_unknowns'], inplace=True)

# now we can drop duplicates and simply keep the first row
# This works because 'False' comes before 'True' in the 'contains_unknowns' column and we want to keep the rows without
# unknown variables if we have information on the same sequence in the same structure in a different row
fasta_ex_df_no_duplicates = fasta_ex_df.drop_duplicates(subset=['gene_name', 'structure_id', 'sequence'], keep='first').copy()

# now we can drop the last column ('contains_unknowns') from both the fasta_ex_df and the
# fasta_ex_df_no_duplicates and sort them according to index!
for df in [fasta_ex_df, fasta_ex_df_no_duplicates]:
    df.drop('contains_unknowns', inplace=True, axis=1)
    df.sort_index(inplace=True)

# now we merge the two dataframes fasta_df and fasta_ex_no_duplicates
combined_df = pd.merge(fasta_df, fasta_ex_df_no_duplicates, how='left', on=['gene_name', 'structure_id', 'sequence'], suffixes=(None, '_ex'))
# we rearrange the columns of the combined dataframe to make the output easily readible
# current cols = ['gene_name', 'structure_id', 'record_id', 'chain_name', 'species',
#                        'description', 'sequence', 'record_id_ex', 'chain_name_ex', 'uniprot_id', 'description_ex']
# we also drop the following columns as we don't need their info in the output: record_id, record_id_ex, chain_name_ex
combined_df = combined_df[['gene_name', 'structure_id', 'chain_name', 'uniprot_id', 'description', 'species', 'description_ex', 'sequence']]

# change back to results dir and write output files with information on all fasta files from all genes
os.chdir(results_dir)

# now, if this is a webrun, we want to read in the data from previous mutafy runs and
# get a slice for all the genes of the current seach (otherwise only new structures will be listed in the output file)
# also, we want to update already existing mutafy data and add new structures to it!
if web_run:
    if exists(f'{mutafy_directory}/04_fasta_info_mutafy.csv'):
        mutafy_fasta_info = pd.read_csv(f'{mutafy_directory}/04_fasta_info_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        # in order to do that, we convert all the values in the fasta_df to strings
        fasta_df = fasta_df.astype('str')
        # now we can concatenate the two dfs
        updated_mutafy_fasta_info = pd.concat([mutafy_fasta_info, fasta_df], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name and then structure id
        updated_mutafy_fasta_info.sort_values(by=['gene_name', 'structure_id', 'record_id'], inplace=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_fasta_info.reset_index(drop=True).to_csv(f'{mutafy_directory}/04_fasta_info_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the resolutions for all the structures for all genes of the current webrun
        fasta_df = updated_mutafy_fasta_info[updated_mutafy_fasta_info.gene_name.isin(list(folder_info.gene_name.values))].reset_index(drop=True)
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        fasta_df.to_csv(f'{mutafy_directory}/04_fasta_info_mutafy.csv', index = False)
        
    if exists(f'{mutafy_directory}/04_fasta_ex_info_mutafy.csv'):
        mutafy_fasta_ex_info = pd.read_csv(f'{mutafy_directory}/04_fasta_ex_info_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        # in order to do that, we convert all the values in the fasta_df to strings
        fasta_ex_df = fasta_ex_df.astype('str')
        # now we can concatenate the two dfs
        updated_mutafy_fasta_ex_info = pd.concat([mutafy_fasta_ex_info, fasta_ex_df], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name and then structure id
        updated_mutafy_fasta_ex_info.sort_values(by=['gene_name', 'structure_id', 'record_id'], inplace=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_fasta_ex_info.reset_index(drop=True).to_csv(f'{mutafy_directory}/04_fasta_ex_info_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the resolutions for all the structures for all genes of the current webrun
        fasta_ex_df = updated_mutafy_fasta_ex_info[updated_mutafy_fasta_ex_info.gene_name.isin(list(folder_info.gene_name.values))].reset_index(drop=True)
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        fasta_ex_df.to_csv(f'{mutafy_directory}/04_fasta_ex_info_mutafy.csv', index = False)
        
    if exists(f'{mutafy_directory}/04_fasta_combined_info_mutafy.csv'):
        mutafy_fasta_combi_info = pd.read_csv(f'{mutafy_directory}/04_fasta_combined_info_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        # in order to do that, we convert all the values in the fasta_df to strings
        combined_df = combined_df.astype('str')
        # now we can concatenate the two dfs
        updated_mutafy_fasta_combi_info = pd.concat([mutafy_fasta_combi_info, combined_df], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name and then structure id
        updated_mutafy_fasta_combi_info.sort_values(by=['gene_name', 'structure_id', 'chain_name'], inplace=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_fasta_combi_info.reset_index(drop=True).to_csv(f'{mutafy_directory}/04_fasta_combined_info_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the resolutions for all the structures for all genes of the current webrun
        combined_df = updated_mutafy_fasta_combi_info[updated_mutafy_fasta_combi_info.gene_name.isin(list(folder_info.gene_name.values))].reset_index(drop=True)
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        combined_df.to_csv(f'{mutafy_directory}/04_fasta_combined_info_mutafy.csv', index = False)

# WRITE RESULTS TO CSV FILES
print(f'>>> writing csv files with information on all genes to {results_dir}') 
print('    >>> writing csv file containing information extracted from fasta files for all genes...')
print('    >>> writing csv file containing information extracted from fasta_ex files for all genes...')
print('    >>> writing csv file containing combined information extracted from fasta and fasta_ex files for all genes...')

fasta_df.to_csv('04_fasta_info.csv', index = False)
fasta_ex_df.to_csv('04_fasta_ex_info.csv', index = False)
combined_df.to_csv('04_fasta_combined_info.csv', index = False)

print('Complete!\n')

# change back to target directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Parsed a total of {fasta_total} fasta files and {fasta_ex_total} fasta_ex files stored across {n_folders} folders.')

print('\nThe following files have been created for each gene and stored in the respective folder:')
print('   o      GENENAME_04_fasta_info.csv                (lists information extracted from all fasta files for this gene/folder)')
print('   o      GENENAME_04_fasta_ex_info.csv           (lists information extracted from all fasta_ex files for this gene/folder)')
print('   o      GENENAME_04_fasta_combined_info.csv (lists combined information extracted from all fasta and fasta_ex files for this gene/folder)')

print('\nThe following files have been created and stored in the Results folder:')
print('   o      04_fasta_info.csv                (lists information extracted from all fasta files for all genes)')
print('   o      04_fasta_ex_info.csv           (lists information extracted from all fasta_ex files for all genes)')
print('   o      04_fasta_combined_info.csv (lists combined information extracted from all fasta and fasta_ex files for all genes)\n\n')

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