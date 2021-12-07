# This script takes two csv files as input, namely:
#      - 01_search_overview_folders.csv
#      - 06_blast_fasta.csv OR 07_blast_two_sequences.csv
# and will
#      - search for pdb files in all folders listed in 01_search_overview_folders.csv
#      - extract info on missing residues / residues which have not been solved in the crystal structure for each pdb structure
#      - combine this information with 06_blast_fasta.csv df
#      - output the following files:
#                - a csv file called 08_unsolved_residues_per_structure.csv listing all unsolved residues in all structures for all genes
#                  (one row for each structure)
#                - a csv files called 08_unsolved_residues_per_chain.csv listing all unsolved residues in all chains of all structures for all genes
#                  (one row for each chain)
#                - a csv file called 08_all_info.csv containing all the information from 06_blast_fasta.csv as
#                  as well as information on unsolved residues extracted from pdb files
# 
#  ----------------------------------------------------------------------------------------------------------------------------------
   
# Set up
import numpy as np
import pandas as pd
import os
from os import listdir
from os.path import isfile, join

from Bio.PDB import *

# ----------------------------------------------------------------------------------------------------------------------------------
target_directory = os.getcwd()
results_dir = f'{target_directory}/Results'

# read in csv files
folders = pd.read_csv(f'{results_dir}/01_search_overview_folders.csv')
# either read in 08_blast_two_sequences.csv or 06_blast_fasta.csv
info = pd.read_csv(f'{results_dir}/07_blast_two_sequences.csv')
# info = pd.read_csv(f'{results_dir}/06_blast_fasta.csv')

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
    # the gene for this row is contained in the folder_name
    # folder_name has the format GENENAME_nstructures
    gene = row.folder_name.split('_')[0]
    
    # change into this folder
    os.chdir(row.full_path)

    # get a list of all pdb files in this folder
    files = [f for f in listdir(row.full_path) if isfile(join(row.full_path, f))]
    pdb_files = [file for file in files if '.pdb' in file]
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
        # we parse the header like so:
        header = parse_pdb_header(pdb)
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
        
        # we also make a df containing one row per chain (insted of one row per structure), which we can later combine more easily with the info df
        for key, value in missing_res_dict.items():
                unsolved_per_chain.loc[len(unsolved_per_chain)] = [gene, pdb[:4], key, value]
        

# change to results directory
os.chdir(results_dir)

# write output to csv files
print('>>> writing csv file containing unsolved residues per structure for all genes...')
print('>>> writing csv file containing unsolved residues per chain for all structures for all genes...\n')
unsolved.to_csv('08_unsolved_residues_per_structure.csv', index=False)
unsolved_per_chain.to_csv('08_unsolved_residues_per_chain.csv', index=False)

        
# to do: combine unsolved df and info df
# idea: loop over column chain_name in info df and find corresponding information for all chains in unsolved_per_chain df

print('>>> combining information on unsolved residues with existing dataframe...')

# we first add an extra columns to the info df which we populate with unsolved residues
# (of all chains in this row/with this sequence in this structure)
info['unsolved_residues_in_structure'] = np.nan

for index, row in info.iterrows():
    # chain_name is in format 'Chain B' or for multiple chains in format 'Chains A, B, C, D, E, F, G, H'
    # we first delete the word Chains and then Chain in case it's only a single chain
    chains = row.chain_name.replace('Chains ', '')
    chains = chains.replace('Chain ', '')
    # now the format of chains is 'B' or A, B, C, D, E, F, G, H' if there are multiple chains
    # in order to get this in a list format, we do
    # for multiple chains:
    if len(chains) > 1:
        chains = chains.split(',')
    # for one chain
    else:
        chains = list(chains)
        
    # now we have a list of all chains with identical sequence in the given structure
    # we now loop over this list and retrieve missing residues from the other df (unsolved_per_chain) for each of them
    # and attach this information to the info df
    # first we create an empty dictionary to populate with all chains (keys) and unsolved residues (values) in this structure
    # which have the same sequence and are thus in one and the same row in the info df
    # (but found across multiple rows in the unsolved_per_chain df)
    unsolved_dict = {}
    for chain in chains:
        # find the unsolved residues for this chain in the unsolved_per_chain df (if there are any)
        # we use a try statement, because if there are no missing residues for this chain, the values will be empty and it will throw an IndexError
        # use structure_id and gene and chain name to identify correct value
        try:
            unsolved_residues = unsolved_per_chain[(unsolved_per_chain.gene == row.gene_name) &
                                               (unsolved_per_chain.structure_id == row.structure_id) &
                                               (unsolved_per_chain.chain == chain)].unsolved_residues_in_chain.values[0]
            # now we append this value and the chain name to the dictionary:
            unsolved_dict[chain] = unsolved_residues
        except IndexError:
            continue

    # now that we've looped over all chains which have the same sequence in this structure, we add the unsolved_dict to the info df
    info.loc[index, 'unsolved_residues_in_structure'] = str(unsolved_dict)

# finally we write the update info df into a csv file called 07_all_info.csv
print('>>> writing csv file called 08_all_info.csv containing all information for all genes...')
info.to_csv('08_all_info.csv', index=False)

# change back to target_directory
os.chdir(target_directory)


print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Parsed a total of {pdb_total} pdb files stored across {len(folders)} folders.')

print('\nThe following files have been created and stored in the Results folder:')
print('   o      08_all_info.csv                                            (lists all information on all structures of all genes (one row = one sequence))')
print('   o      08_unsolved_residues_per_structure.csv     (lists unsolved residues per structure for all genes (one row = one structure))')
print('   o      08_unsolved_residues_per_chain.csv           (lists unsolved residues per chain for all genes (one row = one chain))\n\n')
