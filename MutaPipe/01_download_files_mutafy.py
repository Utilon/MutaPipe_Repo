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

web_run = True # specify if pdb mmcif and fasta files should be stored in separate directory

mutafy_directory = f'{target_directory}/mutafy' # set path to folder where structures will be/are stored
                                            
                                            
# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****    This script takes a csv file containing gene names and corresponding PDB IDs as input and will:
1. create a folder for each gene in the results folder
2. download specified formats (cif, pdb, fasta) into respective folders (atm automatically downloads mmCIF, pdb and FASTA files)
3. outputs a csv file called 01_search_overview_folders listing all the the newly created folders and their contents 
4. outputs a csv file called 01_search_overview_n_structures.csv listing the number of structures retrieved per gene    ***""")

ap.add_argument('-f','--format', nargs='+', required=False, help=f"Specify file format to be downloaded. For mmCif files (.cif) use 'cif' ; for pdb files (.pdb) use 'pdb' ; for fasta files (.fasta) use 'fasta' ; default = {download_format}")
ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in current directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')
ap.add_argument("-w", "--web_run", type=str2bool, required = False, help=f'Indicate whether MutaPipe is run via a webserver (True) or not (False), default = {str(mutafy_directory)}')
ap.add_argument("-m", "--mutafy", required = False, help=f'set path to mutafy directory where information from previous runs is stored, default = {mutafy_directory}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
download_format = download_format if args["format"] == None else args["format"]
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
web_run = web_run if args["web_run"] == None else args["web_run"]
mutafy_directory = mutafy_directory if args["mutafy"] == None else args["mutafy"]

# ----------------------------------------------------------------------------------------------------------------------------------
# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# For the webserver run, some of the data will be in the mutafy_directory:
# create the directory to store downloaded files if it doesn't already exist:
if web_run:
    if not os.path.exists(mutafy_directory):
        os.mkdir(mutafy_directory)
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

# load pdb ids to download
pdb_ids = pd.read_csv(f'{results_dir}/00_search_overview_PDBids.csv', usecols=['gene_name', 'available_structures'], index_col='gene_name')

# if this is a webserver run we also load the pdb ids of all already downloaded and parsed structures from the mutafy directory:
# this won't work if the file doesn't exist yet (first webrun) or if it has been deleted, so we add a try and except statement
if web_run == True:
    try:
        mutafy_pdb_ids = pd.read_csv(f'{mutafy_directory}/00_search_overview_PDBids_mutafy.csv', usecols=['gene_name', 'available_structures'], index_col='gene_name')
    except FileNotFoundError:
        # if the file doesn't exist, we create an empty df so the rest of the code works well
        mutafy_pdb_ids = pd.DataFrame(columns=['gene_name', 'available_structures'])
        mutafy_pdb_ids.set_index('gene_name', inplace=True)

# create an empty list to populate with folder names of all created folders:
created_folders = []

# set variable n_genes: number of genes with available structures to be downloaded
n_genes = len(pdb_ids)

# set variable to count total number of pdb ids for which we download files:
n_structures = 0                  

# create an empty list to populate with gene and number of structures for each gene
list_n_structures = []

# for the webrun, we want to also create a list to keep track of how many *new* structures per gene we have to download
list_new_structures_mutafy = []

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
    
    # create a row with the gene name and the number of available structures for this gene
    # and append row to list_n_structures list
    list_n_structures.append({'gene': gene, 'n_structures': len(found_pdbs)})
    
    # COMPARE WITH DATA FROM PREVIOUS MUTAFY RUNS (IN MUTAFY DATABASE)
    # for the webserver run, we check if the identified PDB ids have already been
    # found/downloaded/parsed in a previous search/run
    if web_run == True:
        # check if the gene is already in the mutafy database (meaning it has been processed in previous webserver runs)
        if gene in mutafy_pdb_ids.index:
            # get a list of pdb ids for this gene which have been processed before by mutafy
            mutafy_pdbs_for_this_gene = ast.literal_eval(mutafy_pdb_ids.loc[gene].available_structures)
            # get a list of only the pdb ids which are new/still have to be processed/downloaded/parsed etc.
            pdb_ids_to_download = [pdb_id for pdb_id in found_pdbs if pdb_id not in mutafy_pdbs_for_this_gene]
            # get length of the list = number of new pdb_ids
            n_new_structures = len(pdb_ids_to_download)
        # if the gene is not in the mutafy database, we still need to define a variable n_new_structures
        # and one pdb_ids_to_download for the script to work
        elif gene not in mutafy_pdb_ids.index:
            # in this case n_new_structures = len(found_pdbs) because all identified PDB ids still have to be downloaded
            n_new_structures = len(found_pdbs)
            pdb_ids_to_download = found_pdbs
        
        # we add these values to the list_new_structures_mutafy which we'll convert to a df
        # and write to a file at the end of the script
        # this file will be used in the following scripts to ensure we only parse new files (and not all files in a folder)
        list_new_structures_mutafy.append({'gene': gene, 'n_structures': len(found_pdbs), 'n_new_structures': n_new_structures, 'new_pdb_ids': pdb_ids_to_download})
        
        # SET UP FOLDER STRUCTURE FOR WEBRUN VS NO WEBRUN 
        # create new folder to save all pdb/mmcif/fasta files found for the query gene:
        # create the folder name for this gene
        folder_name = f'{mutafy_directory}/{gene}_{len(found_pdbs)}structures'
        # get a list of all the currently existing gene folders (if any) to check if a folder for this specific gene has already been created
        existing_folders = [f.path for f in os.scandir(mutafy_directory) if f.is_dir()]
    # for non-webserver-runs the folder name is different, like so:
    else:
        folder_name = f'{results_dir}/{gene}_{len(found_pdbs)}structures'
        # get a list of all the currently existing gene folders (if any) to check if a folder for this specific gene has already been created
        existing_folders = [f.path for f in os.scandir(results_dir) if f.is_dir()]
        
    # it could be the case that a folder for the current gene already exists, but that there are more structures available now
    # in this case we change to the existing folder for the current gene, and rename it to adequately reflect the number of structures available now
    # this ensures that we don't download all already downloaded files again into a new folder
    # we check if any of the existing folders start with the gene name and a '_' (e.g. 'SOD1_')
    # first we get a list of all currently existing folder which fulfill these conditions
    # each element in the list existing_folders is an entire folder path, e.g.
    #  '/scratch/users/k1800109/2023/MutaPipe_Repo/MutaPipe/mutafy/MORC3_6structures'
    # we have to split each folder path at '/' and then check if the last element + one more character correspond to the gene name + '_' (e.g. 'SOD1_')
    if any(f'{gene}_' == existing_folder.split('/')[-1][:len(gene)+1] for existing_folder in existing_folders):
        matching_existing_folders = [existing_folder for existing_folder in existing_folders if existing_folder.split('/')[-1][:len(gene)+1] == f'{gene}_']
        # this list will have at least one element, because we check with any() before we make it and only create it if at least one element in the
        # existing folders list fulfills the criteria
        # if there is 1 element, this is the folder for the current gene and we change to that and rename it
        # if there is more than just 1 entry in the list, we print a warning!!! and then make a new folder / treat it as if there were no matching folders!
        if len(matching_existing_folders) == 1:
            existing_folder_name = matching_existing_folders[0]        
            # if a folder for the current gene exists, we rename it
            os.rename (existing_folder_name, folder_name)
        elif len(matching_existing_folders) > 1:
            print(f'\nWARNING !!! WARNING !!! WARNING!!! \nMultiple folders potentially matching your gene ({gene}) have been identified:\n{[e for e in matching_existing_folders]}\n')
                                  
    # create folder if it doesn't already exist
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    # finally we change to the gene folder (created or renamed) and append the foldername to the list created_folders
    os.chdir(folder_name)
    created_folders.append(folder_name)
     
    
    ## DOWNLOAD THE FILES
    #####################
    # now we download the files - if this is a webserver run, we only download the new ones 
    # First we identify which files already exist in the current folder
    # to do this we get a list of all files
    files = [f for f in listdir(folder_name) if isfile(join(folder_name, f))]
    cif_files = [file for file in files if ('.cif' in file)]
    pdb_files = [file for file in files if ('.pdb' in file)]
    fasta_files = [file for file in files if ('.fasta' in file)]
    
    # for mmCIF files
    if 'cif' in download_format:
        # identify which cif files still need to be downloaded
        if web_run == True:
            cif_to_download = [pdb_id for pdb_id in pdb_ids_to_download if f'{pdb_id}.cif' not in cif_files]
        elif web_run == False:
            cif_to_download = [pdb_id for pdb_id in found_pdbs if f'{pdb_id}.cif' not in cif_files]
        
        if len(cif_to_download) >= 1:
            print(f'\n>>> Initiating download of mmCIF files for {len(cif_to_download)} structures for gene {gene_counter} of {n_genes}: {gene}')
            # we create a PDBList object to download the files with BioPython
            # for mmCIF files
            cifl = PDBList()
            cifl.download_pdb_files(pdb_ids_to_download, file_format='mmCif', pdir=folder_name)
        elif len(cif_to_download) == 0:
            print(f'No new mmCIF structures to be dowloaded for gene {gene_counter} of {n_genes}: {gene}')

    # for PDB files
    if 'pdb' in download_format:        
        # figure out which structures still need to be downloaded and which ones exist
        if web_run == True:
            still_to_be_downloaded = [pdb_id for pdb_id in pdb_ids_to_download if f'{pdb_id}.pdb' not in pdb_files]
        elif web_run == False:
            still_to_be_downloaded = [pdb_id for pdb_id in found_pdbs if f'{pdb_id}.pdb' not in pdb_files]
        
        if len(still_to_be_downloaded) >= 1:
            print(f'\n>>> Initiating download for pdb files for {len(still_to_be_downloaded)} structures for gene {gene_counter} of {n_genes}: {gene}')
            # we create a PDBList object to download the files with BioPython
            pdbl = PDBList()
            pdbl.download_pdb_files(still_to_be_downloaded, file_format='pdb', pdir=folder_name)
        elif len(still_to_be_downloaded) == 0:
            print(f'No new pdb structures to be dowloaded for gene {gene_counter} of {n_genes}: {gene}')
        
        # now we change the names of all pdb files in the folder to be like 4uxy.pdb
        # currently all pdb files have a name like 'pdb4uxy.ent'
        # First, we get an updated list of all pdb files stored in this folder
        files = [f for f in listdir(folder_name) if isfile(join(folder_name, f))]
        new_pdb_files = [f for f in files if '.ent' in f]
        for file in new_pdb_files:
            # the pdb id is the last four letters of the first element of the split
            pdb_id = file.split('.')[0][-4:]
            os.rename(file, pdb_id+'.pdb')

    # for fasta files
    if 'fasta' in download_format:
        # check which fasta files are availalble (if any) and which still need to be downloaded
        # we make a list list with the ones that have to be downloaded
        # if this is a web run, we have to use the list pdb_ids_to_download as a basis (only includes new structures)
        if web_run == True:
            fastas_to_download = [pdb_id for pdb_id in pdb_ids_to_download if f'{pdb_id}.fasta' not in fasta_files]
        # if this isn't a webrun, we have to use the list found_pdbs instead
        elif web_run == False:
            fastas_to_download = [pdb_id for pdb_id in found_pdbs if f'{pdb_id}.fasta' not in fasta_files]
        
        # now we download all fastas which haven't been downloaded already
        if len(fastas_to_download) >= 1:
            print(f'\n>>> Initiating download for fasta files for {len(fastas_to_download)} structures for gene {gene_counter} of {n_genes}: {gene}')
            for pdb_id in fastas_to_download:
                fasta_filename = f'{pdb_id}.fasta'    
                fasta_url = f'https://www.rcsb.org/fasta/entry/{pdb_id.upper()}'
                response = requests.get(fasta_url)
                if response.status_code == 200:
                    print(f'Downloading fasta file for {pdb_id}...')
                    with open(fasta_filename, 'w') as fasta:
                        fasta.write(response.text)
                else:
                    print(response.status_code, response.text)
                    print(f'No fasta file retrieved for {pdb_id}\n')        
        elif len(fastas_to_download) == 0:
            print(f'No new fasta files to be dowloaded for gene {gene_counter} of {n_genes}: {gene}')
            
    print(f'Complete!\n    All corresponding files (format: {download_format}) for {gene} are stored in: \n    {folder_name}\n')
    # remove created folder 'obsolete' (obsolete structure would be stored here if we set obsolete=True for the pdb/cif file download)
    try:
        os.rmdir(folder_name + '/obsolete')
    except:
        pass
    
    # change directory back to results to create next gene folder there
    os.chdir(results_dir) 

# create file with information on folder and contents
# we create an empty list to store rows for a df we create later
list_folders = []
# create lists with filenames of all mmCIF, pdb, and fasta files in this folder
for folder in created_folders:
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    cif_files = [f for f in files if '.cif' in f]
    pdb_files = [f for f in files if '.pdb' in f]
    fasta_files = [f for f in files if '.fasta' in f]
    # sort lists
    cif_files.sort()
    pdb_files.sort()
    fasta_files.sort()
    # we get the gene name and the isolated folder name for the current folder
    # e.g. just "ERBB4_10structures" (i.e. not the full path "/Users/debs/OneDrive - King's College London/2023/Mutafy/MutaPipe_Repo/MutaPipe/mutafy/ERBB4_10structures")
    if web_run == True:
        gene_name = folder.replace((mutafy_directory+'/'), '').split('_')[0] #this is the name of the gene
        isolated_folder_name = folder[len(mutafy_directory)+1:]
    else:
        gene_name = folder.replace((results_dir+'/'), '').split('_')[0] #this is the name of the gene
        isolated_folder_name = folder[len(results_dir)+1:]
    
    # we append data to our list list_folders which we can later convert into a df and write to a file
    # 'folder_name', 'full_path', 'n_mmCIF', 'n_fasta', 'mmCIF_ids', 'fasta_ids'
    list_folders.append({'gene_name': gene_name,
                         'folder_name': isolated_folder_name,
                         'full_path': folder,
                         'n_mmCIF': len(cif_files),
                         'n_pdb': len(pdb_files),
                         'n_fasta': len(fasta_files),
                         'mmCIF_ids': cif_files,
                         'pdb_ids': pdb_files,
                         'fasta_ids': fasta_files})

# now we create dfs from our lists, list_folders and list_n_structures
df_folders = pd.DataFrame(list_folders).astype('str')
df_n_structures   = pd.DataFrame(list_n_structures).astype('str')

# write csv files
df_folders.to_csv('01_search_overview_folders.csv', index = False)
df_n_structures.to_csv('01_search_overview_n_structures.csv', index= False)

# if this is a webrun, we need to convert the list_new_structures_mutafy to a df and write to a file too
# or update the file if it already exists
if web_run:
    df_new_structures_mutafy = pd.DataFrame(list_new_structures_mutafy)
    df_new_structures_mutafy.to_csv(f'{mutafy_directory}/01_new structures_to_be parsed_mutafy.csv', index=False)

# change back to target directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'Complete! \n    downloaded all files (in specified format: {download_format}) for a total of {n_structures} PDB IDs associated with {n_genes} genes\n')
if web_run == True:
    print(f'Folders created/updated in {mutafy_directory}/ :')
    for folder in created_folders:
        print(f'    {folder[len(mutafy_directory)+1:]}')
else:   
    print(f'Folders created/updated in {results_dir}/ :')
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