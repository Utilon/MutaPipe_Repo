# This script takes a csv file (01_search_overview_folders.csv) containing information on the folders where the mmCIF files are stored as input and will:
#      - extract information from each mmCIF file, incl.:
#                - resolution
#                - experimental methods (currently not; maybe add!?)
#                - polypetide sequences (corresponds to chains as shown in PyMOL, e.g. there are missing residues which haven't been solved in the crystal structure)
#                - sequence information in fasta format
#      - outputs the following files:
#                - fasta_ex files for all structures (fasta files created from the mmCIF files)
#                - a csv file called 02_all_resolutions.csv containing the resolutions of all parsed structures for all genes
#                - a csv file called 02_all_poly_seq.csv containing all polypeptide sequences in all structures of all genes
#                - a csv file called 02_structure_info.csv which lists all available header information for each structure for all genes
#  ----------------------------------------------------------------------------------------------------------------------------------
   
# Set up
import pandas as pd
import os
import ast
from os import listdir
from os.path import isfile, join
from os.path import exists
import sys
import argparse
from datetime import datetime

from Bio import SeqIO
from Bio.PDB import *
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder

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

create_search_log = False             # will create a file called search_log.txt with console output if set to True,
                                                     # prints to console if set to False.
target_directory = os.getcwd()      # set target directory (where Results folder is located)
extract_pp = True                         # specify whether to extract the polypeptide sequences or whther to skip this step
delete_files= False                         # specify whether to delete cif files after parsing them or not
web_run = True                             # specify if pdb mmcif and fasta files should be stored in separate directory
mutafy_directory = f'{target_directory}/mutafy' # set path to folder where structures will be/are stored
                                            
                                            
# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****    This script takes a csv file (01_search_overview_folders.csv) containing information on the folders where the mmCIF files are stored as input and will:
1. extract information from each mmCIF file, incl.: resolution, experimental method, polypetide sequences, fasta sequence 
2. outputs the following files: 
(1) fasta_ex files for all structures (fasta files created from the mmCIF files) 
(2) a csv file called GENENAME_02_resolutions.csv per gene containing the resolution for each structure associated with this gene 
(3) a csv file called GENENAME_02_poly_seq.csv per gene containing all the polypeptide sequences in each of the structures associated with this gene 
(4) a csv file called GENENAME_02_structure_info.csv per gene which lists all available header information for each structure of the respective gene 
(5) a csv file called 02_all_resolutions.csv containing the resolutions of all parsed structures for all genes 
(6) a csv file called 02_all_poly_seq.csv containing all polypeptide sequences in all structures of all genes 
(7) a csv file called 02_structure_info.csv which lists all available header information for each structure for all genes     ***""")

ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')
ap.add_argument("-pp", "--polypeptides", type=str2bool, required = False, help=f'Specify whether to extract polypeptide sequence (True) or not (False), default = {str(extract_pp)}')
ap.add_argument("-del", "--delete_files", type=str2bool, required = False, help=f'Specify whether to delete mmCIF files after parsing (True) or not (False), default = {str(delete_files)}')
# ap.add_argument("-s", "--sep_dir", type=str2bool, required = False, help=f'specify if pdb mmcif and fasta files should be stored in separate directory (True) or not (False), default = {str(download_files_to_separate_directory)}')
ap.add_argument("-w", "--web_run", type=str2bool, required = False, help=f'Indicate whether MutaPipe is run via a webserver (True) or not (False), default = {str(mutafy_directory)}')
ap.add_argument("-m", "--mutafy", required = False, help=f'set path to mutafy directory where information from previous runs is stored, default = {mutafy_directory}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
extract_pp = extract_pp if args["polypeptides"] == None else args["polypeptides"]
delete_files = delete_files if args["delete_files"] == None else args["delete_files"]
# download_files_to_separate_directory = download_files_to_separate_directory if args["sep_dir"]   == None else args["sep_dir"]
web_run = web_run if args["web_run"] == None else args["web_run"]
mutafy_directory = mutafy_directory if args["mutafy"] == None else args["mutafy"]
# ----------------------------------------------------------------------------------------------------------------------------------
# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# ----------------------------------------------------------------------------------------------------------------------------------

#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_02.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_02.txt', 'a')

# print nice title
print('===============================================================================')
print('*****    Parsing mmCif Files from the Protein Data Bank for Input Genes    *****')
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
folder_info = pd.read_csv(f'{results_dir}/01_search_overview_folders.csv', usecols=['gene_name', 'folder_name', 'full_path'])

# if this is a web run, we also read in the csv file from script 01 which lists all the new pdb ID's to be downloaded / parsed!
if web_run:
    df_new_structures_to_parse = pd.read_csv(f'{mutafy_directory}/01_new structures_to_be parsed_mutafy.csv')

# define variable n_folders (number of folders to be checked for files) = number of rows in dataframe
n_folders = len(folder_info)

# create variable to keep track of how many cif_files get parsed overall
cif_total = 0

# create list_all_resolutions and list_all_poly_seq and list_all_info to populate with respective information
# and be converted into dfs later on
list_all_resolutions = []
list_all_poly_seq = []
list_all_info = []

# Reading downloaded mmcif files with BioPython
# and extracting resolution and polypeptide sequence from each mmCIF file:
# we loop over the df containing the folder names and the full paths to each folder
# we use this information to parse the files in each folder and extract information
folder_counter = 0
   
for index, row in folder_info.iterrows():
    folder_counter += 1
    gene = row.gene_name
    structure_folder = row.full_path
    # change to the folder containing the files to be parsed
    os.chdir(structure_folder)
    # create list with filenames of all pdb/mmCIF files in this folder
    files = [f for f in listdir(row.full_path) if isfile(join(row.full_path, f))]
    # we get a list of all cif files
    # these are the ones which will be parsed unless it's a webrun and they have been parsed previously 
    cif_files = [f for f in files if '.cif' in f]
    
    # if this is a webrun, we first check if there are any NEW structures to be parsed in this folder/for this gene:
    if web_run:
        # get the list of pdb ids to be parsed        
        new_structures_to_parse = ast.literal_eval(df_new_structures_to_parse[df_new_structures_to_parse.gene == gene].new_pdb_ids.values[0])
        # currently this list contains pdb ids, but we need mmCif filenames (pdb_id.cif)
        # so we do the following:
        cif_files = [f'{pdb_id}.cif' for pdb_id in new_structures_to_parse]
    
    # sort list
    cif_files.sort()
    # if there are no new structures to be parsed, we can continue to the next gene/folder
    if len(cif_files) == 0:
        print(f'\nNo new mmCif files to be parsed for {gene} (gene {folder_counter} of {n_folders})')
        continue
    
    print(f"""\nStarting to parse identified mmCif files for {gene} (gene {folder_counter} of {n_folders}):                   {len(cif_files)} mmCIF files""")
        
    # add number of cif files in folder to total number of cif files:
    cif_total += len(cif_files)
            
    # for each mmCIF file to be imported,
    # set the filename and structure_id that BioParser needs to import mmCIF files and load structure objects       
    cif_counter = 0
    for cif_file in cif_files:
        cif_counter += 1        
        structure_id = cif_file.replace('.cif', '')          # get the PDB id for each file to be parsed
        
        # we can also extract a FASTA file directly from a mmCIF with the sequence in FASTA format:
        print(f'        >>> creating FASTA file extracted from mmCIF file for {structure_id}')
        seq = SeqIO.parse(cif_file, 'cif-seqres')
        SeqIO.write(seq, f'{structure_id}_ex.fasta', 'fasta')
        
        # To load structures for cif files, we first create an MMCIFParser object:
        parser = MMCIFParser(QUIET=True)
        # load structure object:
        print(f'    Getting {gene} structure object for {structure_id}                                                                ({cif_counter} of {len(cif_files)} from mmCIF files for {gene})')
        # added try and except statement to capture cases where the cif file doesn't exist
        # (no clue why this should be the case but we will print a WARNING if so)
        try:
            structure = parser.get_structure(structure_id, cif_file)
        except FileNotFoundError:
            print(f'WARNING! No file {cif_file} exists')
            continue
        
        # extract header information:
        # =====================
        # extract information from created BioPython structure object
        print(f'        >>> extracting header information mmCIF file for {structure_id}')        
        resolution = structure.header["resolution"]
        structure_method = structure.header['structure_method']
        deposition_date = structure.header['deposition_date']
        structure_name = structure.header['name']
        classification = structure.header['head']
        
        # not all structures have resolutions (depending on structure method)!
        # as we want to sort our output later on according to resolution, we are going to
        # replace all missing values with 999
        if not resolution:
            resolution = 999

        # create a new entry in both the list_all_info,  and list_all_resolutions for this structure
        # add new row to each list (will later be converted to df)
        list_all_info.append({'gene': gene,
                             'structure_id': structure_id,
                             'resolution': resolution,
                             'structure_method': structure_method,
                             'deposition_date': deposition_date,
                             'structure_name': structure_name,
                             'classification' : classification})
        
        list_all_resolutions.append({'gene': gene,
                                     'structure_id': structure_id,
                                     'resolution': resolution})

        # get polypeptide sequences for all polypeptides in current structure
        # ==================================================
        # The polypeptide sequences correspond to the sequence as seen in pyMOL, i.e. with gaps/missing residues
        # we only want to do this if the argument --polypeptides is set to True, therefore I added an if statement
        if extract_pp == True:
            # to temporarily store polypeptide sequences:
            # create an empty dictionary for polypeptide sequences:
            all_poly_seqs = {}
                
            print(f'        >>> getting polypeptides sequences from mmCIF file for {structure_id}')
            PolypeptideBuilder = PPBuilder()
            polypeptides = PolypeptideBuilder.build_peptides(structure)
            # Sometimes the polypeptides cannot be extraced properly with the PolypeptideBuilder,
            # (don't know why, but this happens with all the KIF5A structures, in this case polypeptides is an empty list)
            if not polypeptides:
                print(' Could not build polypeptide sequences with PPBuilder.')
                
            # create new dictionary entry in all_poly_seqs dictionary in the format {pdb id:[seq1, seq2...]}
            # create entry (empty list to be populated) for key=structure_id
            all_poly_seqs[structure_id] = []
                
            # add all polypeptides to list in all_poly_seqs[structure_id]
            # if polypeptides is empty because the PPBuilder did not work properly, this code will not do anything as there are no elements in the list to loop over.
            for pp in polypeptides:
                seq = pp.get_sequence()
                # add sequence to all_poly_seqs dictionary (important to add as a string):
                all_poly_seqs[structure_id].append(str(seq))
                
            # Append data to the list_all_poly_seq
            list_all_poly_seq.append({'gene': gene,
                                  'structure_id': structure_id,
                                  'n_polypeptides': len(all_poly_seqs[structure_id]),
                                  'polypeptides': all_poly_seqs[structure_id]})
                        
        # now that we have finished processing the structure object extracted from the cif file,
        # we can delete this cif file to save space on the disk
        if delete_files == True:
            os.remove(cif_file)    
    
    print(f'Complete!\n    All mmCIF files for {gene} have been parsed!')
    
##### CONVERT TO DFS AND WRITE CSV FILES
# Convert the lists to dfs so we can concatencate with existing mutafy data and write them to csv files
df_all_resolutions = pd.DataFrame(list_all_resolutions)
df_all_info = pd.DataFrame(list_all_info)
df_all_poly_seq = pd.DataFrame(list_all_poly_seq)
# convert list columns to string, otherwise concatenation later doesn't work
df_all_poly_seq = df_all_poly_seq.astype({'polypeptides':'str'})

# now, if this is a webrun, we want to read in the data from previous mutafy runs and
# get a slice for all the genes of the current seach (otherwise only new structures will be listed in the output file)
# also, we want to update already existing mutafy data and add new structures to it!
if web_run:
    if exists(f'{mutafy_directory}/02_all_resolutions_mutafy.csv'):
        mutafy_resolutions = pd.read_csv(f'{mutafy_directory}/02_all_resolutions_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        updated_mutafy_resolutions = pd.concat([mutafy_resolutions, df_all_resolutions], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name and then structure id
        updated_mutafy_resolutions.sort_values(by=['gene', 'structure_id'], inplace=True, ignore_index=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_resolutions.to_csv(f'{mutafy_directory}/02_all_resolutions_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the resolutions for all the structures for all genes of the current webrun
        df_all_resolutions = updated_mutafy_resolutions[updated_mutafy_resolutions.gene.isin(list(folder_info.gene_name.values))]
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        df_all_resolutions.to_csv(f'{mutafy_directory}/02_all_resolutions_mutafy.csv', index = False)

    if exists(f'{mutafy_directory}/02_structure_info_mutafy.csv'):
        mutafy_structure_info = pd.read_csv(f'{mutafy_directory}/02_structure_info_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        updated_mutafy_structure_info = pd.concat([mutafy_structure_info, df_all_info], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name and then structure id
        updated_mutafy_structure_info.sort_values(by=['gene', 'structure_id'], inplace=True, ignore_index=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_structure_info.to_csv(f'{mutafy_directory}/02_structure_info_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the resolutions for all the structures for all genes of the current webrun
        df_all_info = updated_mutafy_structure_info[updated_mutafy_structure_info.gene.isin(list(folder_info.gene_name.values))]
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        df_all_info.to_csv(f'{mutafy_directory}/02_structure_info_mutafy.csv', index = False)
        
# update mutafy polypeptide outputs from all runs
    if exists(f'{mutafy_directory}/02_all_poly_seq_mutafy.csv'):
        mutafy_pp_seq = pd.read_csv(f'{mutafy_directory}/02_all_poly_seq_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        updated_mutafy_pp_seq = pd.concat([mutafy_pp_seq, df_all_poly_seq], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name and then structure id
        updated_mutafy_pp_seq.sort_values(by=['gene', 'structure_id'], inplace=True, ignore_index=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_pp_seq.to_csv(f'{mutafy_directory}/02_all_poly_seq_mutafy.csv', index=False)
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the resolutions for all the structures for all genes of the current webrun
        df_all_poly_seq = updated_mutafy_pp_seq[updated_mutafy_pp_seq.gene.isin(list(folder_info.gene_name.values))]
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        df_all_poly_seq.to_csv(f'{mutafy_directory}/02_all_poly_seq_mutafy.csv', index = False)
    
# change back to Results directory
os.chdir(results_dir)

# write resolutions of all structures for all genes to csv file
df_all_resolutions.to_csv('02_all_resolutions.csv', index = False)

# write all parsed info of all structures for all genes to csv file
df_all_info.to_csv('02_structure_info.csv', index = False)

# write polypeptide sequences of all structures for all genes to csv file:
if extract_pp == True:
    df_all_poly_seq.to_csv('02_all_poly_seq.csv', index = False)

# change back to target directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Succesfully parsed a total of {cif_total} mmCIF files stored across {n_folders} folders.')

print('\nThe following files have been created for each gene and stored in the respective folder:')
if extract_pp == True:
    print('   o      GENENAME_02_poly_seq.csv          (lists all polypeptide sequences for each structure of the respective gene)')
print('   o      GENENAME_02_resolutions.csv        (lists the resolution for each structure of the respective gene; 999 indicates a missing value)')
print('   o      GENENAME_02_structure_info.csv   (lists all available header information for each structure of the respective gene)')

print('\nThe following files have been created and stored in the Results folder:')
if extract_pp == True:
    print('   o      02_all_poly_seq.csv                (lists all polypeptide sequences for each structure of all genes)')
print('   o      02_all_resolutions.csv            (lists the resolution for each structure of all genes; 999 indicates a missing value)')
print('   o      02_structure_info.csv            (lists all available header information for each structure for all genes)\n\n')

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