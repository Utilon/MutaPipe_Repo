# This script takes a csv file containing gene names and corresponding PDB IDs as input and will:
#      - create a folder for each gene in the results folder
#      - download specified formats (cif, pdb, fasta) into respective folders 
#      - outputs a csv file called 01_search_overview_folders listing all the the newly created folders and their contents
    
#  ----------------------------------------------------------------------------------------------------------------------------------
   
# Set up
import pandas as pd
import ast
import requests
import os
from os import listdir
from os.path import isfile, join
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
# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory
# ----------------------------------------------------------------------------------------------------------------------------------

#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_01.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_01.txt', 'a')
    
# print nice title
print('===============================================================================')
print('*****    Downloading Files from the Protein Data Bank for Input Genes    *****')
print('===============================================================================\n')

# print script name to console/log file
print(f'script name: {script_name}')

# store current date and time in an object and print to console / write to log file
start_time = datetime.now()
print(f'start: {start_time}\n')
# ----------------------------------------------------------------------------------------------------------------------------------


# load pdb ids to download
results_dir = f'{target_directory}/Results'
pdb_ids = pd.read_csv(f'{results_dir}/00_search_overview_PDBids.csv', usecols=['gene_name', 'available_structures'], index_col='gene_name')

# create an empty list to populate with folder names of all created folders:
created_folders = []

# set variable n_genes: number of genes with available structures to be downloaded
n_genes = len(pdb_ids)

# set variable to count total number of pdb ids for which we download files:
n_structures = 0

# create empty df to populate with gene and number of structures:
df_n_structures = pd.DataFrame(columns=['gene', 'n_structuctures'])

# Initiate for loop to iterate over pdb_ids df and download all mmCIF, pdb and fasta files for each gene
gene_counter = 0

for gene, structures in pdb_ids.iterrows():
    gene_counter += 1    
    # gene is the gene name and the row index of the pdb_ids df
    # structures is a Pandas series.
    # In order to get a list of all the structures in the right format:
    for structure in structures:                                   # there is only one entry per structures, it's a list representation in string format
        found_pdbs = ast.literal_eval(structure)           # use ast.literal_eval to convert the string into a list
    
    # update total number of structures
    n_structures += len(found_pdbs)
    
    # append a row with the gene name and the number of available structures for this gene to the df_n_structures df
    df_n_structures.loc[len(df_n_structures)] = [gene, len(found_pdbs)]

    # create new folder to save all pdb/mmcif/fasta files found for the query gene:
    folder_name = f'{results_dir}/{gene}_{len(found_pdbs)}structures'
#     folder_name_add_on = 1 #number will be added to folder name if folder name already exists # I commented this because, we actually don't want to make a new folder if the folder already exists, I think that's better

#         commented the while loop because we don't make a new folder, instead I added an if statement here
#     while os.path.exists(folder_name):
#         folder_name = f'{results_dir}/{gene}_{len(found_pdbs)}structures_{folder_name_add_on}'
#         folder_name_add_on += 1
    if os.path.exists(folder_name):
        os.chdir(folder_name)
        created_folders.append(folder_name)
    else:
         os.makedirs(folder_name)
         os.chdir(folder_name)
         created_folders.append(folder_name)
            
    # download all mmCIF files to newly created folder
    # if-statement added, so mmCif files are only downloaded if specified (default)
    if 'cif' in download_format:
        print(f'\n>>> Initiating download of mmCIF files for {len(found_pdbs)} structures for gene {gene_counter} of {n_genes}: {gene}')
        # we create a PDBList object to download the files with BioPython
        # for mmCIF files
        cifl = PDBList()
        cifl.download_pdb_files(found_pdbs, file_format='mmCif', pdir=folder_name)
        
    # for PDB files
    # if-statement added, so pdb files are only downloaded if specified (default)
    if 'pdb' in download_format:
        print(f'\n>>> Initiating download for pdb files for {len(found_pdbs)} structures for gene {gene_counter} of {n_genes}: {gene}')
        # first we check which pdb files are available (if any) in this folder, so we don't download the same files again:
        # create list with filenames of all pdb files in this folder
        files = [f for f in listdir(folder_name) if isfile(join(folder_name, f))]
        pdb_files = [file for file in files if ('.pdb' in file)]
        # in order to check if we have downloaded the same structures again, we need to reconvert the name
        # of these files to their original name, which is in the format 'pdb4uxy.ent'
        # (they are currently in the format ''4uxy.pdb')
        for file in pdb_files:
            orig_name = 'pdb' + file[:4] + '.ent'
            os.rename(file, orig_name)
        
        # now we can download the structures and will automatically
        # get a warning if the structure is already available in our directory as we do for the mmCif files
        # (this doesn't work for the fasta files as I don't use BioPython for it, so I implemented it manually, see further below)
        # we create a PDBList object to download the files with BioPython        
        pdbl = PDBList()
        pdbl.download_pdb_files(found_pdbs, file_format='pdb', pdir=folder_name)
    
        # currently all pdb files have a name like 'pdb4uxy.ent'
        # so in order to rename the files to pdb_id.pdb (e.g. 4uxy.pdb), we do the following:
        # First, we get a list of all pdb files stored in this folder
        files = [f for f in listdir(folder_name) if isfile(join(folder_name, f))]
        pdb_files = [f for f in files if '.ent' in f]
        for file in pdb_files:
            # the pdb id is the last four letters of the first element of the split
            pdb_id = file.split('.')[0][-4:]
            os.rename(file, pdb_id+'.pdb')
    
    # now we download the fasta file for this structure
    # if-statement added, so fasta files are only downloaded if specified (default)
    if 'fasta' in download_format:
        print(f'\n>>> Initiating download for fasta files for {len(found_pdbs)} structures for gene {gene_counter} of {n_genes}: {gene}')
        # first we check which fasta files are available (if any) in this folder, so we don't download the same files again:
        # create list with filenames of all fasta files in this folder
        files = [f for f in listdir(folder_name) if isfile(join(folder_name, f))]
        fasta_files = [file for file in files if ('.fasta' in file)]
        
        # we loop over the pdb ids and check if each fasta file exists or not
        for pdb_id in found_pdbs:
            fasta_filename = f'{pdb_id}.fasta'
            if fasta_filename in fasta_files:
                print(f'Fasta file exists: {os.getcwd()}/{fasta_filename}')
                continue
            else:                                    
                fasta_url = f'https://www.rcsb.org/fasta/entry/{pdb_id.upper()}'
                response = requests.get(fasta_url)
                if response.status_code == 200:
                    print(f'Downloading fasta file for {pdb_id}...')
                    with open(fasta_filename, 'w') as fasta:
                        fasta.write(response.text)
                else:
                    # If there is no data, print status code and response
                    print(response.status_code, response.text)
                    print(f'No fasta file retrieved for {pdb_id}\n')        
                
    print(f'Complete!\n    All corresponding files (format: {download_format}) for {gene} are stored in: \n    {folder_name}\n')
    # remove created folder 'obsolete' (obsolete structure would be stored here if we set obsolete=True for the pdb/cif file download)
    try:
        os.rmdir(folder_name + '/obsolete')
    except:
        continue
    # change directory back to results to create next gene folder there
    os.chdir(results_dir) 

    
# change back to results directory
os.chdir(results_dir)    

# create df with information on folder and contents
df_folders = pd.DataFrame(columns=['gene_name', 'folder_name', 'full_path', 'n_mmCIF', 'n_pdb', 'n_fasta', 'mmCIF_ids', 'pdb_ids', 'fasta_ids'])

# create lists with filenames of all mmCIF and fasta files in this folder
for folder in created_folders:
    gene_name = folder.replace((results_dir+'/'), '').split('_')[0] #this is the name of the gene
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    cif_files = [f for f in files if '.cif' in f]
    pdb_files = [f for f in files if '.pdb' in f]
    fasta_files = [f for f in files if '.fasta' in f]
    # sort lists
    cif_files.sort()
    pdb_files.sort()
    fasta_files.sort()
    # append data to df_folders df:
    # 'folder_name', 'full_path', 'n_mmCIF', 'n_fasta', 'mmCIF_ids', 'fasta_ids'
    df_folders.loc[len(df_folders)] = [gene_name, folder[len(results_dir)+1:], folder, len(cif_files), len(pdb_files), len(fasta_files), cif_files, pdb_files, fasta_files]
    
# write csv files
df_folders.to_csv('01_search_overview_folders.csv', index = False)

df_n_structures.to_csv('01_search_overview_n_structures.csv', index= False)

# change back to target directory
os.chdir(target_directory)


print('\n============================== Summary ================================================\n')
print(f'Complete! \n    downloaded all files (in specified format: {download_format}) for a total of {n_structures} PDB IDs associated with {n_genes} genes\n')
print(f'Folders created/updated in {results_dir}/ :')
folder_counter = 1
n_folders = len(created_folders)
for folder in created_folders:
    print(f'    {folder[len(results_dir)+1:]}')

print('\nThe following files have been created:')
print('   o      01_search_overview_folders.csv              (lists all the the created/updated folders and their contents)')
print('   o      01_search_overview_n_structures.csv          (lists number of structures retrieved per gene)\n\n')


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