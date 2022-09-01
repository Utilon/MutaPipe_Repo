# This script takes a csv file (01_search_overview_folders.csv) containing information on the folders where the mmCIF files are stored as input and will:
#      - extract information from each mmCIF file, incl.:
#                - resolution
#                - experimental methods (currently not; maybe add!?)
#                - polypetide sequences (corresponds to chains as shown in PyMOL, e.g. there are missing residues which haven't been solved in the crystal structure)
#                - sequence information in fasta format
#      - outputs the following files:
#                - fasta_ex files for all structures (fasta files created from the mmCIF files)
#                - a csv file called GENENAME_02_resolutions.csv per gene containing the resolution for each structure associated with this gene
#                - a csv file called GENENAME_02_poly_seq.csv per gene containing all the polypeptide sequences in each of the structures associated with this gene.
#                - a csv file called GENENAME_02_structure_info.csv per gene which lists all available header information for each structure of the respective gene
#                - a csv file called 02_all_resolutions.csv containing the resolutions of all parsed structures for all genes
#                - a csv file called 02_all_poly_seq.csv containing all polypeptide sequences in all structures of all genes
#                - a csv file called 02_structure_info.csv which lists all available header information for each structure for all genes
#  ----------------------------------------------------------------------------------------------------------------------------------
   
# Set up
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
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

create_search_log = False     # will create a file called search_log.txt with console output if set to True,
                                            # prints to console if set to False.
target_directory = os.getcwd()    # set target directory (where Results folder is located)
extract_pp = True                 # specify whether to extract the polypeptide sequences or whther to skip this step
delete_files= False                  # specify whether to delete cif files after parsing them or not
download_files_to_separate_directory = True # specify if pdb mmcif and fasta files should be stored in separate directory
                                            
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
ap.add_argument("-s", "--sep_dir", type=str2bool, required = False, help=f'specify if pdb mmcif and fasta files should be stored in separate directory (True) or not (False), default = {str(download_files_to_separate_directory)}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
extract_pp = extract_pp if args["polypeptides"] == None else args["polypeptides"]
delete_files = delete_files if args["delete_files"] == None else args["delete_files"]
download_files_to_separate_directory = download_files_to_separate_directory if args["sep_dir"]   == None else args["sep_dir"]

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

# Read in data from csv file
folder_info = pd.read_csv(f'{results_dir}/01_search_overview_folders.csv', usecols=['gene_name', 'folder_name', 'full_path'])

# define variable n_folders (number of folders to be checked for files) = number of rows in dataframe
n_folders = len(folder_info)

# create variable to keep track of how many cif_files get parsed overall
cif_total = 0

# create df_all_resolutions and df_all_poly_seq and df_all_info to populate with respective information
df_all_resolutions = pd.DataFrame(columns=['gene', 'structure_id', 'resolution'])
df_all_poly_seq = pd.DataFrame(columns=['gene', 'structure_id'])
df_all_info = pd.DataFrame(columns=['gene', 'structure_id', 'resolution', 'structure_method', 'deposition_date', 'structure_name', 'classification'])

# Reading downloaded mmcif files with BioPython
# and extracting resolution and polypeptide sequence from each mmCIF file:
# we loop over the df containing the folder names and the full paths to each folder
# we use this information to parse the files in each folder and extract information
folder_counter = 0

for index, row in folder_info.iterrows():
    folder_counter += 1
    gene = row.gene_name      
    
    # change to the folder containing the files to be parsed
    os.chdir(row.full_path)      
        
    # create list with filenames of all pdb/mmCIF files in this folder
    files = [f for f in listdir(row.full_path) if isfile(join(row.full_path, f))]
    cif_files = [f for f in files if '.cif' in f]
    # sort list
    cif_files.sort()
        
    print(f"""\nChecking folder content (folder {folder_counter} of {n_folders}) for gene {gene}:                   {len(cif_files)} mmCIF files""")
    
    # add number of cif files in folder to total number of cif files:
    cif_total += len(cif_files)
    
    # to temporarily store polypeptide sequences:
    # create an empty dictionary for polypeptide sequences:
    all_poly_seqs = {}
        
    # for each mmCIF file to be imported,
    # set the filename and structure_id that BioParser needs to import mmCIF files and load structure objects       
    cif_counter = 0
    for cif_file in cif_files:
        cif_counter += 1        
        structure_id = cif_file.replace('.cif', '')          # get the PDB id for each file to be parsed
        # To load structures for cif files, we first create an MMCIFParser object:
        parser = MMCIFParser(QUIET=True)
        
        # load structure object:
        print(f'    Getting {gene} structure object for {structure_id}                                                                ({cif_counter} of {len(cif_files)} from mmCIF files for {gene})')
        structure = parser.get_structure(structure_id, cif_file)
                
        # we can also extract a FASTA file with the sequence in FASTA format:
        print(f'        >>> creating FASTA file extracted from mmCIF file for {structure_id}')
        seq = SeqIO.parse(cif_file, 'cif-seqres')
        SeqIO.write(seq, f'{structure_id}_ex.fasta', 'fasta')
        
        # now that we have extracted the structure object as well as the fasta file from the cif file,
        # we can delete this cif file to save space on the disk
        if delete_files == True:
            os.remove(cif_file)      
        
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
        
        # create a new entry in the df_all_info,  and df_all_resolutions for this structure
        # add new row to each df
        df_all_info.loc[len(df_all_info)] = [gene, structure_id, resolution, structure_method, deposition_date, structure_name, classification]
        df_all_resolutions.loc[len(df_all_resolutions)] = [gene, structure_id, resolution]


        # get polypeptide sequences for all polypeptides in current structure
        # ==================================================
        # The polypeptide sequences correspond to the sequence as seen in pyMOL, i.e. with gaps/missing residues
        # we only want to do this if the argument --polypeptides is set to True, therefore I added an if statement
        if extract_pp == True:
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
            counter=0
            for pp in polypeptides:
                counter += 1
                seq = pp.get_sequence()
                # add sequence to all_poly_seqs dictionary:
                all_poly_seqs[structure_id].append(seq)
        
    print(f'Complete!\n    All mmCIF files for {gene} have been parsed!')
                
    # create csv files for each gene with resolution and all info of all structures for that gene
    print(f'        >>> saving csv files containing resolution and all parsed info for all {gene} structures')
    # if we use a separate directory for the downloaded files, we have to write them to the results forlder instead, so they do
    # not get overwritten
    if download_files_to_separate_directory == True:
        df_all_resolutions[df_all_resolutions.gene == gene].to_csv(f'{results_dir}/{gene}_02_resolutions.csv', index = False)
        df_all_info[df_all_info.gene == gene].to_csv(f'{results_dir}/{gene}_02_structure_info.csv', index = False)
    else:
        # if the option is false, we store the csv files in the current directory (where the structure files are also stored)
        df_all_resolutions[df_all_resolutions.gene == gene].to_csv(f'{gene}_02_resolutions.csv', index = False)
        df_all_info[df_all_info.gene == gene].to_csv(f'{gene}_02_structure_info.csv', index = False)

    # now, only if --polypeptides is set to true, we want to
    # create csv file for each gene with sequences of all polypeptides of all structures for that gene
    if extract_pp == True:
        # first create a pandas dataframe from dictionary containing all sequences for all structures of the current gene
        # first we create a list of column names for our dataframe:
        # create column names for each polypeptide: pp1 - ppn
        if polypeptides:
            column_names = []
            # find longest list in dictionary to create the same amount of columns in pd dataframe:
            max_num_of_polyseqs_per_structure = 0
            for value in all_poly_seqs.values():
                if len(value) > max_num_of_polyseqs_per_structure:
                    max_num_of_polyseqs_per_structure = len(value)
                    
            counter = 0
            for i in range(max_num_of_polyseqs_per_structure):
                counter += 1
                col_name = 'pp'+str(counter)
                column_names.append(col_name)
            # create dataframe
            df = pd.DataFrame.from_dict(all_poly_seqs, orient='index', columns=column_names).reset_index()
            df = df.rename(columns={'index':'structure_id'})
            # sort dataframe by 'id' (alphabetically)
            df = df.sort_values(by=['structure_id'])
            # save dataframe to csv file
            print(f'        >>> saving csv file containing all polypeptide sequences for all {gene} structures')
            if download_files_to_separate_directory == True:
                df.to_csv(f'{results_dir}/{gene}_02_poly_seq.csv', index = False)
            else:
                df.to_csv(f'{gene}_02_poly_seq.csv', index = False)
            
            # to add it to the overall df df_all_poly_seq, we first add a column with the gene name
            df['gene'] = gene
            # now we can merge the two dfs with an outer merge.
            # we cannot simply append the df because they have different numbers of columns (due to different numbers of polypeptides)
            df_all_poly_seq = df_all_poly_seq.merge(df, how='outer')

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