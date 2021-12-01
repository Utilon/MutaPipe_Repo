# This script takes a csv file containing gene names and corresponding PDB IDs as input and will:
#      - create a folder for each gene in the results folder
#      - download specified formats (cif, pdb, fasta) into respective folders 
#        (atm automatically downloads mmCIF, pdb and FASTA files --> update later with argparse to be able to specify desired file format via terminal)
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

# load pdb ids to download
target_directory = os.getcwd()
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

# Initiate for loop to iterate over pdb_ids df and download all mmCIF and fasta files for each gene
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
    folder_name_add_on = 1 #number will be added to folder name if folder name already exists
        
    while os.path.exists(folder_name):
        folder_name = f'{results_dir}/{gene}_{len(found_pdbs)}structures_{folder_name_add_on}'
        folder_name_add_on += 1
    else:
         os.makedirs(folder_name)
         os.chdir(folder_name)
         created_folders.append(folder_name)
            
    # download all mmCIF files to newly created folder
    print(f'\n>>> downloading mmCIF, pdb and fasta files for {len(found_pdbs)} structures for gene {gene_counter} of {n_genes}: {gene}...')
    # we create a PDBList object to download the files with BioPython
    # for mmCIF files
    cifl = PDBList()
    cifl.download_pdb_files(found_pdbs, file_format='mmCif', pdir=folder_name)
    
    # for PDB files
    pdbl = PDBList()
    pdbl.download_pdb_files(found_pdbs, file_format='pdb', pdir=folder_name)
    
    # currently all downloaded pdb files have a name like 'pdb4uxy.ent'
    # so in order to rename the files to pdb_id.pdb (e.g. 4uxy.pdb), we do the following:
    # First, we get a list of all pdb files stored in this folder
    files = [f for f in listdir(folder_name) if isfile(join(folder_name, f))]
    pdb_files = [f for f in files if '.ent' in f]
    for file in pdb_files:
        # the pdb id is the last four letters of the first element of the split
        pdb_id = file.split('.')[0][-4:]
        os.rename(file, pdb_id+'.pdb')
    
    # now we download the fasta file for this structure    
    for pdb_id in found_pdbs:
        fasta_url = f'https://www.rcsb.org/fasta/entry/{pdb_id.upper()}'
        response = requests.get(fasta_url)
        if response.status_code == 200:
            print(f'Downloading fasta file for {pdb_id}...')
            with open(f'{pdb_id}.fasta', 'w') as fasta:
                fasta.write(response.text)
        else:
            # If there is no data, print status code and response
            print(response.status_code, response.text)
            print(f'No fasta file retrieved for {pdb_id}\n')        
        
    print(f'Complete!\n    All mmCIF, pdb and fasta files for {gene} are stored in: \n    {folder_name}')
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
print(f'Complete! \n    downloaded all mmCIF, pdb and fasta files for a total of {n_structures} PDB IDs associated with {n_genes} genes\n')
print(f'New folders created in {results_dir}/ :')
folder_counter = 1
n_folders = len(created_folders)
for folder in created_folders:
    print(f'    {folder[len(results_dir)+1:]}')

print('\nThe following files have been created:')
print('   o      01_search_overview_folders.csv              (lists all the the newly created folders and their contents)')
print('   o      01_search_overview_n_structures.csv          (lists number of structures retrieved per gene)\n\n')