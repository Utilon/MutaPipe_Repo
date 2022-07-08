# Script to get best n structures per mutation
# This script takes the following csv files as input:
#       - 02_structure_info.csv
#       - 03_unsolved_residues_per_chain.csv
#       - 05_blastp_results.csv
#       - 06_b_ClinVar_Annotations.csv
# it will perform the following:
#      - combine the information in the 3 dfs (according to PDBid and chain)
#      - filter out sequences which are shorter than a given percentage of the reference sequence (set variable relative_sequence_length)
#      - filter out sequences whose best hsp covers less than a given percentage of the reference sequence (set variable hsp_coverage)
#      - sort/filter the df in order to get:
#               - n best structures (best resolution) for all single amino acid variants (SAVs) (structures with only this one mutation and no other mutations)
#               - n best structures (best resolution) for all unique combinations of mutations available in the pdb
#               - n best structures (best resolution) for any specific mutation, regardless of other mutations in the same structure
#               - all wildtype structures (defined as HSP covering 99% of reference sequence, 100% similarity, no mismatches)
#      - add all availalable ClinVar annotations to all three n_best_structure tables/dfs
#      - add additional info (validation report) from the PDB API to the all output tables/dfs 
# and outputs the following files:
#      - In each respective gene folder:
#               - GENENAME_07_best_structures_per_SAV.csv
#                 lists n best structures for each SAV (one mutation per structure) for this gene (incl. ClinVar annotations)
#               - GENENAME_07_best_structures_all_unique_combinations.csv
#                 lists n best structure for all unique mismatch combinations for this gene (incl. ClinVar annotations)
#               - GENENAME_07_best_structures_any_mutation.csv
#                 lists n best structures for any mismatch regardless of other mismatches in the same structure (incl. ClinVar annotations)
#               - GENENAME_07_wildtype_structures.csv
#                 lists all available WT structures for this gene
#      - In the results folder:
#               - 07_best_structures_per_SAV.csv
#                  lists best structure for each point mutation (one mutation per structure) in all genes (incl. ClinVar annotations)
#               - 07_best_structures_all_unique_combinations.csv
#                  lists best structure for all unique mismatch combinations for all genes (incl. ClinVar annotations)
#               - 07_best_structures_any_mutation.csv
#                  lists best structure for any mismatch regardless of other mismatches in all genes (incl. ClinVar annotations)
#               - 07_wildtype_structures.csv
#                 lists all available WT structures for this gene

# ===========================================================================================

# Set up
import pandas as pd
import numpy as np
import os
import ast
import sys
import argparse
import requests
from datetime import datetime

# get this script's name:
script_name = os.path.basename(__file__)

# suppress pandas SettingWithCopyWarning
pd.options.mode.chained_assignment = None  # default='warn'

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

# define a function to make number inputs between 0.1 and 1.0 possible for hsp_coverage and relative_sequence_length
def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x
    
# set default values for arguments we want to implement
# we have to do this here if we want to print the default values in the help message

create_search_log = False     # will create a file called search_log.txt with console output if set to True,
                                            # prints to console if set to False.
relative_sequence_length = 0.5               # filter out sequences which are shorter than a given percentage of
                                                              # the reference sequence (set variable relative_sequence_length)
                                                              # e.g. 50% (0.5) or 10% (0.1) of the reference sequence

hsp_coverage = 0.1                                # filter out sequences whose best hsp covers less than a given
                                                              # percentage of the reference sequence
                                                              #  e.g. 50% (0.5) or 10% (0.1) of the reference sequence
                                                              
target_directory = os.getcwd()    # set target directory (where Results folder is located)

n_best_structures = 5 # number of top structures you want to include in the output table (1 for the best, 2 for the two best available one) per sequence/mismatch

exclude_unsolved_mismatches = False # indicate if structures where the mismatch of interest is not solved in the crystal structure should be excluded (True) or not (False)
                                                                                    
# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****     Script to get best n structures per mutation. This script takes the following csv files as input:
- 02_structure_info.csv
- 03_unsolved_residues_per_chain.csv
- 05_blastp_results.csv
- 06_b_ClinVar_Annotations.csv
It will perform the following:
(1) combine the information in the 3 dfs (according to PDBid and chain)
(2) filter out sequences which are shorter than a given percentage of the reference sequence (set variable relative_sequence_length)
(3) filter out sequences whose best hsp covers less than a given percentage of the reference sequence (set variable hsp_coverage)
(4) sort/filter the df in order to get:
    (a) n best structures (best resolution) for all single amino acid variants (SAVs) (structures with only this one mutation and no other mutations)
    (b) n best structures (best resolution) for all unique combinations of mutations available in the pdb
    (c) n best structures (best resolution) for any specific mutation, regardless of other mutations in the same structure
    (d) all wildtype structures (defined as HSP covering 99% of reference sequence, 100% similarity, no mismatches)
(5) add all availalable ClinVar annotations to all three n_best_structure tables/dfs
(6) add additional info (validation report) from the PDB API to the all output tables/dfs 
(7) output the following files:
- In each respective gene folder:
(a) GENENAME_07_best_structures_per_SAV.csv (lists n best structures for each SAV (one mutation per structure) for this gene (incl. ClinVar annotations))
(b) GENENAME_07_best_structures_all_unique_combinations.csv (lists n best structure for all unique mismatch combinations for this gene (incl. ClinVar annotations))
(c) GENENAME_07_best_structures_any_mutation.csv (lists n best structures for any mismatch regardless of other mismatches in the same structure (incl. ClinVar annotations))
(d) GENENAME_07_wildtype_structures.csv (lists all available WT structures for this gene)
- In the results folder:
(a) 07_best_structures_per_SAV.csv (lists best structure for each point mutation (one mutation per structure) in all genes (incl. ClinVar annotations))
(b) 07_best_structures_all_unique_combinations.csv (lists best structure for all unique mismatch combinations for all genes (incl. ClinVar annotations))
(c) 07_best_structures_any_mutation.csv (lists best structure for any mismatch regardless of other mismatches in all genes (incl. ClinVar annotations))
(d) 07_wildtype_structures.csv (lists all available WT structures for this gene)    ***""")

ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')
ap.add_argument("-rsl", "--relative_sequence_length", type=restricted_float, required = False, help=f'filter out sequences shorter than a given percentage of the reference sequence, default = {str(relative_sequence_length)}')
ap.add_argument("-cov", "--hsp_coverage", type=restricted_float, required = False, help=f'filter out sequences whose best hsp covers less than a given percentage of the reference sequence, default = {str(hsp_coverage)}')
ap.add_argument("-n_best", "--n_best_structures", type=int , required = False, help=f'Specify number of top structures per sequence/variant to be included in final output, default = {str(n_best_structures)}')
ap.add_argument("-e", "--exclude_unsolved_mismatches", type=str2bool, required = False, help=f'indicate whether to exclude cases where the mismatch of interest is not solved in the crystal structure (True) or not (False), default = {str(exclude_unsolved_mismatches)}')
args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
relative_sequence_length = relative_sequence_length if args["relative_sequence_length"] == None else args["relative_sequence_length"]
hsp_coverage  = hsp_coverage if args["hsp_coverage"]   == None else args["hsp_coverage"]
n_best_structures = n_best_structures if args["n_best_structures"] == None else args["n_best_structures"]
exclude_unsolved_mismatches = exclude_unsolved_mismatches if args["exclude_unsolved_mismatches"] == None else args["exclude_unsolved_mismatches"]

# ----------------------------------------------------------------------------------------------------------------------------------
# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# ----------------------------------------------------------------------------------------------------------------------------------
#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_07.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_07.txt', 'a')
    
# print nice title
print('===============================================================================')
print('*****    Get the Best Structures for Each Unique Sequence in the PDB for Input Genes    *****')
print('===============================================================================\n')

# print script name to console/log file
print(f'script name: {script_name}')

# store current date and time in an object and print to console / write to log file
start_time = datetime.now()
print(f'start: {start_time}\n')

# ----------------------------------------------------------------------------------------------------------------------------------
# Define functions
def add_empty_cols(df):
    new_df = df.reindex(columns=([col for col in df.columns] +
                                 ['accession', 'title', 'variant_type', 'protein_change', 'aliases', 'clinical significance',
                                  'last_evaluated', 'review_status', 'associated_traits', 'dbs_and_accessions']), fill_value='')
    return new_df

# Define a function that can change one letter AA codes to three letter AA codes and vice versa
def change_aa_code(one_or_three_letter_code):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA':'X'}
      
    if len(one_or_three_letter_code) == 3:
        updated_code = d[one_or_three_letter_code]
    elif len(one_or_three_letter_code) == 1:
        updated_code = list(d.keys())[list(d.values()).index(one_or_three_letter_code)]
    else:
        print('Not a valid one-/three-letter AA code, no conversion applied')
        updated_code = one_or_three_letter_code
    return updated_code

# Define a function to add clinvar data to dfs
def add_clinvar_annotations(df):
    
    # reset index of df! (drop old index with drop=True)
    df.reset_index(drop=True, inplace=True)
    
    for index, row in df.iterrows():
        gene = row.gene_name
        # get a slice of the clinvar_annotations df with variants for only the current gene
        clinvar_annotations_this_gene = clinvar_annotations[clinvar_annotations.gene == gene]
        
        # we add a try and except statement to use the mismatch_of_interest column if available
        # and otherwise we work with the columns mismatch_substitutions and close_mismatch_substitutions
        # (only the case in the unique_combi df)
        try:
            mismatch = row.mismatch_of_interest
            # now we make two lists with all potential mismatches and aliases from clinvar from this df slice for our current mismatch
            potential_corresponding_clinvar_mismatches = [x for x in clinvar_annotations_this_gene.protein_change if mismatch in str(x)]
            potential_corresponding_clinvar_aliases = [x for x in clinvar_annotations_this_gene.aliases if mismatch in str(x)]
        except AttributeError:
            all_mismatches = ast.literal_eval(row.mismatch_substitutions) + ast.literal_eval(row.close_mismatch_substitutions)
            potential_corresponding_clinvar_mismatches = [x for x in clinvar_annotations_this_gene.protein_change if any(mismatch in str(x) for mismatch in all_mismatches)]
            potential_corresponding_clinvar_aliases = [x for x in clinvar_annotations_this_gene.aliases if any(mismatch in str(x) for mismatch in all_mismatches)]

        # we can use these lists to get a slice of the clinvar_annotations df with all potential correspoonding clinvar mismatches
        potential_mismatches_df_slice = clinvar_annotations_this_gene[clinvar_annotations_this_gene.protein_change.isin(potential_corresponding_clinvar_mismatches)]
        potential_aliases_df_slice = clinvar_annotations_this_gene[clinvar_annotations_this_gene.aliases.isin(potential_corresponding_clinvar_aliases)]

        # we concatenate the two df slices to get a df with all the clinvar annotations for this variant/the current mismatch
        corresponding_clinvar_annotations = pd.concat([potential_mismatches_df_slice, potential_aliases_df_slice])
        
        # now we can use the information from the concatenated df and append it to the current row in the df (the one we are looping over)
        # ideally there is only one row in the df, i.e. only one variant/alias with this name
        # (or if there is no row at all, we don't have to append anything)
        # if len(corresponding_clinvar_annotations) > 1, we need to figure out another way to add the information to the current row
        if len(corresponding_clinvar_annotations) == 1:
            # we fill the rows 'accession' to 'dbs_and_accessions' of the df with the the values from all but the first two columns of the one entry in the corresponding_clinvar_annotations df
            df.loc[index, 'accession':'dbs_and_accessions'] = corresponding_clinvar_annotations.iloc[0, 2:]                
        elif len(corresponding_clinvar_annotations) > 1:
            # if there is more than one row, we want to get the entire series (each column) from corresponding_clinvar_annotations and write that into the cells in the df
            # sometimes this doesn't work, so we add a try and except statement
            # (the values we want to write into the cells don't come as lists but as ndarrays, then we get an error saying:
            # '"Must have equal len keys and value " ValueError: Must have equal len keys and value when setting with an ndarray'
            # I think I can avoid this problem if we take use str in front of every variable
            # so I have removed the try and except statement again
            df.loc[index, 'accession':'dbs_and_accessions'] = [str(corresponding_clinvar_annotations.accession.values),
                                                                                str(corresponding_clinvar_annotations.title.values),
                                                                                str(corresponding_clinvar_annotations.variant_type.values),
                                                                                str(corresponding_clinvar_annotations.protein_change.values),
                                                                                str(corresponding_clinvar_annotations.aliases.values),
                                                                                str(corresponding_clinvar_annotations['clinical significance'].values),
                                                                                str(corresponding_clinvar_annotations.last_evaluated.values),
                                                                                str(corresponding_clinvar_annotations.review_status.values),
                                                                                str(corresponding_clinvar_annotations.associated_traits.values),
                                                                                str(corresponding_clinvar_annotations.dbs_and_accessions.values)]
    return df

# define function for GET request to PDB API
def get_data(url):
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        # If there is no data, print status code and response
        print(response.status_code, response.text)
        print(f'No data file retrieved for input url {url}\n')
        return None

# define a function which takes a pdb id as input and generates
# a csv file with all the relevant information from the PDB API
def get_validation_report(pdb_id):
    # construct urls to query PDB API
    # see https://data.rcsb.org/#rest-api
    # define base url
    base_url = 'https://data.rcsb.org/rest/v1/core'

    # define all specific query urls
    entry_url = base_url + f'/entry/{pdb_id.upper()}'
    # we need to specify the assembly number as well, we take the first one
    assembly = '1'
    assembly_url = base_url + f'/assembly/{pdb_id.upper()}/{assembly}'
    
    # use the GET function and the urls to get data from the PDB API
    entry_data = get_data(entry_url)
    assembly_data = get_data(assembly_url)

    # EXTRACTING INFORMATION FROM THE ENTRY DATA
    # ======================================
    # all attributes listed here:
    # https://data.rcsb.org/data-attributes.html#entry-attributeshttps://data.rcsb.org/data-attributes.html#entry-attributes

    # I think pdbx_vrpt_summary is the validation report for protein structures
    # get data as shown in experimental snapshot on pdb
    # some of these keys might not be available for all structures (e.g. Cryo-EM structures or NMR structures)
    # thus we add a try and exept statement for each key
    try:
        resolution = entry_data['pdbx_vrpt_summary']['pdbresolution']
    except KeyError:
        resolution = np.NaN
    try:
        r_free = entry_data['pdbx_vrpt_summary']['pdbrfree']                                           
    except KeyError:
        r_free = np.NaN
    try:    
        r_work = entry_data['pdbx_vrpt_summary']['pdb_r']
    except KeyError:
        r_work = np.NaN
    try:
        r_observed = entry_data['pdbx_vrpt_summary']['dcc_r']
    except KeyError:
        r_observed = np.NaN
    try:
        clashscore = entry_data['pdbx_vrpt_summary']['clashscore']
    except KeyError:
        clashscore = np.NaN
    try:
        ramachandran_outliers = entry_data['pdbx_vrpt_summary']['percent_ramachandran_outliers']
    except KeyError:
        ramachandran_outliers = np.NaN
    try:
        sidechain_outliers = entry_data['pdbx_vrpt_summary']['percent_rotamer_outliers']
    except KeyError:
        sidechain_outliers = np.NaN
    try:
        rsrz_outliers = entry_data['pdbx_vrpt_summary']['percent_rsrzoutliers']
    except KeyError:
        rsrz_outliers = np.NaN

    # EXTRACTING INFORMATION FROM THE ASSEMBLY DATA
    # ======================================
    # all attributes listed here:
    # https://data.rcsb.org/data-attributes.html#entry-attributeshttps://data.rcsb.org/data-attributes.html#entry-attributes
    assembly_info = assembly_data['rcsb_assembly_info']
    # to get
    polymer_composition = assembly_data['rcsb_assembly_info']['polymer_composition']

    structure_symmetry = assembly_data['rcsb_struct_symmetry'][0]
    # to get
    symmetry = structure_symmetry['type']
    stoichiometry = structure_symmetry['stoichiometry']
    oligomeric_state = structure_symmetry['oligomeric_state']
    
    # ----------------------------
    # include a link to pdf file of the validation report, this has many great figures which I cannot
    # include in the output tables, and is easier to read overall anyway
    # not sure that works for all structures, but for the ones I checked the link seems to be constructed as follows:
    pdf_full_val_report = f'https://files.rcsb.org/pub/pdb/validation_reports/{pdb_id[1:3]}/{pdb_id}/{pdb_id}_full_validation.pdf'

    # Create df to populate with extracted variables
    structure_info = pd.DataFrame({'structure_id' : [pdb_id], 
                                           'resolution' : [resolution],
                                           'r_free' : [r_free],
                                           'r_work' : [r_work],
                                           'r_observed' : [r_observed],
                                           'clashscore' : [clashscore],
                                           'ramachandran_outliers' : [ramachandran_outliers],
                                           'sidechain_outliers' : [sidechain_outliers],
                                           'RSRZ_outliers' : [rsrz_outliers],
                                           'polymer_composition' : [polymer_composition],
                                           'oligomeric_state' : [oligomeric_state],
                                           'symmetry' : [symmetry],
                                           'stoichometry' : [stoichiometry],
                                           'full_validation_report' : [pdf_full_val_report]})
    
    return structure_info



# ----------------------------------------------------------------------------------------------------------------------------------
#  read in data
# ===========
#       - 02_structure_info.csv
#       - 03_unsolved_residues_per_chain.csv
#       - 05_blastp_results.csv
#       - 06_b_ClinVar_Annotations.csv
structure_info = pd.read_csv(f'{results_dir}/02_structure_info.csv')
unsolved_per_chain = pd.read_csv(f'{results_dir}/03_unsolved_residues_per_chain.csv')
blastp_results = pd.read_csv(f'{results_dir}/05_blastp_results.csv')
clinvar_annotations = pd.read_csv(f'{results_dir}/06_b_ClinVar_Annotations.csv')

#  combine dfs
# ===========
# first, we combine the blastp_results df with the unsolved_per_chain df as follows:
# idea: loop over column chain_name in blastp_results df and find corresponding information for all chains in unsolved_per_chain df
# we first add an extra columns to the blastp_results df which we populate with unsolved residues
# (of all chains in this row/with this sequence in this structure)
blastp_results['unsolved_residues_in_structure'] = np.nan
# we also add columns for number and percent of unsolved residues
blastp_results['n_unsolved_residues'] = np.nan
blastp_results['percent_unsolved_residues'] = np.nan

for index, row in blastp_results.iterrows():
    # chain_name is in format 'Chain B' or for multiple chains in format 'Chains A, B, C, D, E, F, G, H'
    # we first delete the word Chains and then Chain in case it's only a single chain
    chains = row.chain_name.replace('Chains ', '')
    chains = chains.replace('Chain ', '')
    # now the format of chains is 'B' or 'A, B, C, D, E, F, G, H' if there are multiple chains
    # in order to get this in a list format, we do
    # for multiple chains:
    if len(chains) > 1:
        chains = chains.split(',')
    # for one chain
    else:
        chains = list(chains)
    # finally we make sure there are no spaces in the chain names,
    # e.g. Chains A, B will be a list like ['A', ' B'] and we need to remove the space in ' B'
    chains = [chain.strip() for chain in chains]
        
    # now we have a list of all chains with identical sequence in the given structure
    # we now loop over this list and retrieve missing residues from the other df (unsolved_per_chain) for each of them
    # and attach this information to the blastp_results df
    # first we create an empty dictionary to populate with all chains (keys) and unsolved residues (values) in this structure
    # which have the same sequence and are thus in one and the same row in the blastp_results df
    # (but found across multiple rows in the unsolved_per_chain df)
    unsolved_dict = {}
    n_unsolved_res_per_chain = {}
    percent_unsolved_per_chain = {}
    for chain in chains:
        # find the unsolved residues for this chain in the unsolved_per_chain df (if there are any)
        # we use a try statement, because if there are no missing residues for this chain, the values will be empty and it will throw an IndexError
        # use structure_id and gene and chain name to identify correct value
        try:
            unsolved_residues = unsolved_per_chain[(unsolved_per_chain.gene == row.gene_name) &
                                               (unsolved_per_chain.structure_id == row.structure_id) &
                                               (unsolved_per_chain.chain == chain)].unsolved_residues_in_chain.values[0]
            # currently this is a list in string format, so we convert it to a list:
            unsolved_residues = ast.literal_eval(unsolved_residues)
            # now we append this value and the chain name to the dictionary:
            unsolved_dict[chain] = unsolved_residues
            # we also get the length of the list unsolved_residues which corresponds to the number of unsolved res in the chain
            # we add this to the n_unsolved_res_per_chain dict
            n_unsolved_res_per_chain[chain] = len(unsolved_residues)
            # we also get the percentage of unsolved residues in the chain (relative to the sequence length)
            percent_unsolved_per_chain[chain] = round(len(unsolved_residues) / len(row.sequence), 3)
        except IndexError:
            # we append an empty list instead if there is no data on unsolved residues for this chain
            unsolved_dict[chain] = []
            n_unsolved_res_per_chain[chain] = np.nan
            percent_unsolved_per_chain[chain] = np.nan
            continue

    # now that we've looped over all chains which have the same sequence in this structure, we add the unsolved_dict to the blastp_results df
    blastp_results.loc[index, 'unsolved_residues_in_structure'] = str(unsolved_dict)
    # we also add the number and the percentage of unsolved residues to the blastp_results df
    blastp_results.loc[index, 'n_unsolved_residues'] = str(n_unsolved_res_per_chain)
    blastp_results.loc[index, 'percent_unsolved_residues'] = str(percent_unsolved_per_chain)

# now the blastp_results df contains all corresponding unsolved residues and more information on them (percentage /number)!
# However, the blast_p_results df does not yet include the resolutions of the structures, so we get the resolutions from
# a previous output (structure_info) by combining the two dataframes:
structure_info.rename(columns={'gene': 'gene_name'}, inplace=True)
df = blastp_results.merge(structure_info, how='left', on=['gene_name', 'structure_id'])    


# FILTER RELATIVE SEQUENCE LENGTH
# filter out structures which are shorter a given percentage of the reference sequence
# --> relative_sequence_length, e.g. 50% (0.5) or 10% (0.1) of the reference sequence
# ==================================================================================
print(f'Relative sequence length threshold set to {relative_sequence_length}:')
print(f'        only strucures with a sequence length of at least {relative_sequence_length*100}%  of the reference sequence will be included in output\n')
df = df[df.sequence.apply(lambda x: len(x)) > relative_sequence_length]

# FILTER RELATIVE HSP COVERAGE
# filter out structures which cover less than hsp_coverage, e.g. 50% (0.5) or 10% (0.1) of the reference sequence
# ==================================================================================
print(f'Hsp coverage threshold set to {hsp_coverage}:')
print(f'        only strucures with an hsp covering at least {hsp_coverage*100}%  of the reference sequence will be included in output\n')
df = df[df.hsp_length / df.alignment_length > hsp_coverage]
# the table lists the best hsp per structure (in case there are multiple hsps for this strucure's alignment with the reference structure)
# we want to create a filter to exclude very short sequences which only cover a small part of the reference sequence
# the alignment_length in the df is always the same as the length of the reference sequence
# so by dividing the hsp_length / alignment_length, we get the a percentage for the coverage of this hsp for the reference sequence
# (e.g. hsp covers 75 % of reference sequence)
# with the filter we can exclude all hsps which are below a certain coverage
# we set this threshold to 10% for now
# note: for certain proteins, e.g. FUS, the reference sequence is very long (FUS: 526 AAs) and the hsps of sequences in the structures are much shorter
# for these genes a lower coverage will have to be set if no structures get identified otherwise
# (this applies in all cases where it's difficult to crystallise the entire protein and hence only small fragments are available, I would assume)

# FILTER WILDTYPE
# ==================================================================================
# The WT should have
#      - 100% similarity to the canonical sequence,
#      - no mismatches compared to the canonical sequence,
#      - the best hsp should cover 99% of the canonical sequence
#        (ideally 100%, but keep it at 99% in case of changes in conventions such as position numbering etc., e.g. like in the case of SOD1)
df_wt = df[(df.similarity == 1.0) & (df.mismatches_incl_gaps == 0) & (df.hsp_length/df.alignment_length >= 0.99)]

# we order the df_wt first by gene name, then by hsp length, then by similarity, then by resolution
df_wt = df_wt.sort_values(by=['gene_name', 'hsp_length', 'similarity', 'resolution'], ascending=[True, False, False, True])

# we keep all the available wilttype structures in the df (don't drop any) 
# write final df to file at the end of the script (see below)


# GET N BEST STRUCTURES FOR ALL SAVS
#  =================================================================================
# next we create a df that contains the n best structure for all single amino acid variants (SAVs; only one mismatch/close_mismatch)
single_mutation = df[df.close_mismatches + df.mismatches_excl_gaps == 1]

# we order the single mutation df first by gene name, then by normal mismatches, by close mismatches, and finally by resolution
single_mut_sorted = single_mutation.sort_values(by=['gene_name', 'mismatch_substitutions', 'close_mismatch_substitutions', 'resolution'])

# we reset the index of our df single_mut_sorted
single_mut_sorted.reset_index(drop = True, inplace = True)

# Filter rows
# in order to get the aa index for each mutation so we are able to sort the df accordingly, we do the following:
# make mismatch of interest col in df single_mut_sorted 
# first we make a column called mismiatch_of_interest
# single_mut_sorted.loc[:,'mismatch_of_interest'] = None # this line does not work with KL 24.6.2022
# try this instead
single_mut_sorted['mismatch_of_interest'] = None  # this works with KL 24.6.2022
# we also make a column called mismatch_solved_in_crystal_structure where we indicate if the mismatch
# of interest has been solved in the structure/chain (True) or not (False)
single_mut_sorted['mismatch_solved_in_crystal_structure'] = None
# to fill the column with values, we loop over the df (couldn't figure out how else to do it)
for index, row in single_mut_sorted.iterrows():
    close_mis = ast.literal_eval(row.close_mismatch_substitutions)
    mis = ast.literal_eval(row.mismatch_substitutions)
    all_mismatches = mis +close_mis
    single_mut_sorted.loc[index, 'mismatch_of_interest'] = all_mismatches[0]
    
# now we use mismatch_of_interest column to create extra col in single_mut_sorted
# containing the AA index of each mismatch (position)
single_mut_sorted.loc[:,'aa_index'] = single_mut_sorted.mismatch_of_interest.apply(lambda x: int(x[1:-1]))

# Now we check if the mismatch of interest is solved in the crystal structure or not:
# we get the mismatch of interest and check if it's in the unsolved residues
# for this, we get the unsolved residues and change the amino acid code from 3 to 1 letter code
for index, row in single_mut_sorted.iterrows():
    this_mismatch = row.mismatch_of_interest
    # get a dict with all unsolved residues per chain
    # in format {'A': ['LYS283', 'ARG184', ...], 'B' : [..., ...]}
    unsolved_residues = ast.literal_eval(row.unsolved_residues_in_structure)
    # create an empty dict to populate with transformed values, i.e. changed AA codes
    one_letter_unsolved_residues = {}
    for key, value in unsolved_residues.items():
        one_letter_unsolved_residues[key] = [change_aa_code(res[:3])+res[3:] for res in value]
    # now we can check if our mismatch_of_interest is listed in the one_letter_unsolved_residues
    # mismatch of interest in format 'T162A', thus we take all but the last character of the string
    # we want to create a dictionary with keys = chains and a Boolean value indicating if the mismatch of interest
    # is solved in the crystal structure (True) or not (False)
    solved_or_not ={}
    for key, value in one_letter_unsolved_residues.items():
        if this_mismatch[:-1] in value:
            # if this_mismatch is in the unsolved residue list, the mismatch of interest has NOT been solved in
            # the crystal structure, so we output False
            solved_or_not[key] = False
        else:
            # otherwise the mismatch of interest has been solved in the crystal structure,
            # and we output True
            solved_or_not[key] = True
        # now we can add the solved_or_not dict to the df in the column 'mismatch_solved_in_crystal_structure'
        single_mut_sorted.loc[index, 'mismatch_solved_in_crystal_structure'] = str(solved_or_not)
        
# Add option to be able to exclude structures where the mismatch of interest is unsolved in crystal structure
# (missing atomic coordinates)     
if exclude_unsolved_mismatches == True:
    # we want to drop all rows where the mismatch of interest is NOT solved in the crystal structure
    # if it's solved in any chain in the corresponding row, the word 'True' will be contained
    # somewhere in the column mismatch_solved_in_crystal_structure'
    for index, row in single_mut_sorted.iterrows():
        # we drop rows if the mismatch has NOT been solved, i.e. True is not in the column mismatch_solved_in_crystal_structure
        if ' True' not in row.mismatch_solved_in_crystal_structure:
            single_mut_sorted.drop(labels=index, axis=0, inplace = True)
    # now that we dropped all the rows for which the mismatch of interest has not been solved in the crystal structure,
    # we can reset the index of our df again
    single_mut_sorted.reset_index(drop = True, inplace = True)
    

# FILTER n_best_structures
# we want to get the first n rows for each gene and each mismatch of interest
# For this, we create an empty df which we populate with the correct data in the loop below:
n_best_structures_single_mut = pd.DataFrame()
        
# get a slice of the single_mut_sorted df for only one gene at a time:
# loop over unique gene_names in single_mut_sorted:
for gene in single_mut_sorted.gene_name.unique():
    # get a slice containing only data of this gene
    single_mut_gene_slice = single_mut_sorted[single_mut_sorted.gene_name == gene]
    # now we can use this slice to get a smaller slice of each mismatch
    for mismatch in single_mut_sorted.mismatch_of_interest.unique():
        single_mut_mismatch_slice = single_mut_gene_slice[single_mut_gene_slice.mismatch_of_interest == mismatch]
        # now sort this slice according to resolution and take the first n rows 
        n_best_single_mut_for_this_gene = single_mut_mismatch_slice.sort_values(by='resolution')[:n_best_structures]
        # and append it to the df we created just before this loop (n_best_structures_single_mut)
        n_best_structures_single_mut = n_best_structures_single_mut.append(n_best_single_mut_for_this_gene)

# now we can sort the df according to the aa_index
# the following line throws an Error (KeyError) if the df is empty (if there are no SAVs)
# e.g. in the case of KL (24.6.2022)
# we thus we add a try and except statement
try:
    df_SAV = n_best_structures_single_mut.sort_values(by=['gene_name', 'aa_index', 'mismatch_of_interest', 'resolution'])
except KeyError:
    # if this throws an Error it's because the df is empty, so we don't sort the df, and use the other df with all relevant column names
    df_SAV = single_mut_sorted
# reset index of df
df_SAV.reset_index(drop = True, inplace = True)
# write final output to file at the end of script (with other outputs)

# GET N BEST STRUCTURES FOR ANY UNIQUE COMBINATION OF MUTATIONS
# =======================================================================================
# get the best structure for any unique combination of mutations per gene
# in order to do that we first create a new column in the df which contains a list of all the mismatch substitutions
# (close mismatches and normal mismatches)
# we make two lists of the mismatch and close_mismatch column and add them to each other
# we also check if the variable is a string. missing values are floats, so they get ignored as we write an empty list instead
df['all_mismatches'] = df.mismatch_substitutions.apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else []) + df.close_mismatch_substitutions.apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else [])

# create a list of all unique gene names, so we can loop over a slice of the df containing only structures of this gene
# we do this because the mismatches for one gene are comparable, but not for multiple genes (e.g. A5E is not the same mutation if it is in two differenct proteins and cannot be summed up in one row)
unique_genes = df.gene_name.unique()

# create empty df to populate with best structure per unique combination:
df_unique_combi = pd.DataFrame(columns=df.columns)

# create empty df to populate with best structure per unique mutation (regardless of other mutations)
best_structure_any_mutation = pd.DataFrame(columns=df.columns)
# this df also needs an extra column called mismatch_of_interest
best_structure_any_mutation['mismatch_of_interest'] = None
# and a column called mismatch_solved_in_crystal_structure where we indicate if the mismatch
# of interest has been solved in the structure/chain (True) or not (False)
best_structure_any_mutation['mismatch_solved_in_crystal_structure'] = None

# now we loop over the unique genes and take a slice of the df containing only these rows (with this gene)
for gene in unique_genes:
    print(f'Filtering data to get the best structures for gene {gene}')
    df_gene = df[df.gene_name == gene].copy()

    # now we use the df_gene for the rest of the loop
    # we create a list of all unique mismatch and close mismatch combinations using the new column called all_mismatches
    # in order to use the unique function, we have to convert the column all_mismatches back into a string format
    # (lists are not hashable!)
    df_gene['all_mismatches'] = df_gene.all_mismatches.apply(lambda x: str(x))
    unique_mismatch_combinations = [element for element in list(df_gene.all_mismatches.unique())]

    # now that we have all unique mismatch combinations, we can filter our df
    # to get only the structures with this particular combination
    # we can then filter again in order to get the n best structures with the highest resolution for each combination!
    # initiate for loop
    for combination in unique_mismatch_combinations:
        # we extract a slice of the df only containing the rows with this mismatch combination
        df_this_combi = df_gene[df_gene.all_mismatches == combination]
        
        # sort by resolution
        df_this_combi_sorted = df_this_combi.sort_values(by='resolution')
        
        # FILTER n_best_structures
        # we want to get the first n rows for each df_this_combi_sorted
        best_n_structures_this_combi = df_this_combi_sorted[:n_best_structures]
        
        # we now append the best_n_structures_this_combi to the df_unique_combi df
        df_unique_combi = df_unique_combi.append(best_n_structures_this_combi)
    # reset index of our df
    df_unique_combi.reset_index(drop = True, inplace = True)
        
    # GET N BEST STRUCTURE FOR ANY MUTATION (still in the loop)
    # =========================================================================================
    # get the best structure for any mutation regardless of whether there are other mutations (also per gene, still in the loop)

    # reconvert the all_mismatches column in df_gene from string to list
    df_gene['all_mismatches'] = df_gene.all_mismatches.apply(lambda x: ast.literal_eval(x))
    
    # we get a list of all unique mismatches (not mismatch combinations, but single residue mutations) for the current gene
    unique_mismatches = []
    
    # we loop over the all_mismatches columns,
    # each row contains a list with all mismatches in this structure
    for mismatches in df_gene.all_mismatches:
        # we loop over all the mismatches in this structure
        for mismatch in mismatches:
            # we append this mismatch to unique_mismatches if it's not already there
            if mismatch not in unique_mismatches:
                unique_mismatches.append(mismatch)

    # we have to convert the all_mismatches column back to a string (it's currently a list)
    df_gene['all_mismatches'] = df_gene.all_mismatches.apply(lambda x: str(x))

    # now we loop over the list of unique mismatches for this gene and get the best structure for each of them
    for mismatch in unique_mismatches:
        # first we get a df of all structures which contain this mismatch (regardless of other mismatches in the structure)
        this_mismatch = df_gene[df_gene.all_mismatches.str.contains(mismatch)]
        # now we sort the this_mismatch df to get the best structures (highest resolution)
        # sort by  resolution only
        this_mismatch_sorted = this_mismatch.sort_values(by=['resolution'])
        this_mismatch_sorted.reset_index(drop = True, inplace = True)
        
        # we add information on this mismatch to the new column we created earlier
        this_mismatch_sorted['mismatch_of_interest'] = mismatch
        
        # Now we check if the mismatch of interest is solved in the crystal structure or not:
        # we get the mismatch of interest and check if it's in the unsolved residues
        # for this, we get the unsolved residues and change the amino acid code from 3 to 1 letter code
        for index, row in this_mismatch_sorted.iterrows():
            this_mismatch = row.mismatch_of_interest
            # get a dict with all unsolved residues per chain
            # in format {'A': ['LYS283', 'ARG184', ...], 'B' : [..., ...]}
            unsolved_residues = ast.literal_eval(row.unsolved_residues_in_structure)
            # create an empty dict to populate with transformed values, i.e. changed AA codes
            one_letter_unsolved_residues = {}
            for key, value in unsolved_residues.items():
                one_letter_unsolved_residues[key] = [change_aa_code(res[:3])+res[3:] for res in value]
            # now we can check if our mismatch_of_interest is listed in the one_letter_unsolved_residues
            # mismatch of interest in format 'T162A', thus we take all but the last character of the string
            # we want to create a dictionary with keys = chains and a Boolean value indicating if the mismatch of interest
            # is solved in the crystal structure (True) or not (False)
            solved_or_not ={}
            for key, value in one_letter_unsolved_residues.items():
                if this_mismatch[:-1] in value:
                    # if this_mismatch is in the unsolved residue list, the mismatch of interest has NOT been solved in
                    # the crystal structure, so we output False
                    solved_or_not[key] = False
                else:
                    # otherwise the mismatch of interest has been solved in the crystal structure,
                    # and we output True
                    solved_or_not[key] = True
                # now we can add the solved_or_not dict to the df in the column 'mismatch_solved_in_crystal_structure'
                this_mismatch_sorted.loc[index, 'mismatch_solved_in_crystal_structure'] = str(solved_or_not)
                
        # Add option to be able to exclude structures where the mismatch of interest is unsolved in crystal structure
        # (missing atomic coordinates)     
        if exclude_unsolved_mismatches == True:
            # we want to drop all rows where the mismatch of interest is NOT solved in the crystal structure
            # if it's solved in any chain in the corresponding row, the word 'True' will be contained
            # somewhere in the column mismatch_solved_in_crystal_structure'
            for index, row in this_mismatch_sorted.iterrows():
                # we drop rows if the mismatch has NOT been solved, i.e. True is not in the column mismatch_solved_in_crystal_structure
                if ' True' not in row.mismatch_solved_in_crystal_structure:
                    this_mismatch_sorted.drop(labels=index, axis=0, inplace = True)
            # now that we dropped all the rows for which the mismatch of interest has not been solved in the crystal structure,
            # we can reset the index of our df again
            this_mismatch_sorted.reset_index(drop = True, inplace = True)      
        
        # FILTER n_best_structures        
        # we want to get the first n rows of the this_mismatch_sorted df
        best_n_structures_this_mismatch = this_mismatch_sorted[:n_best_structures]        
        
        # we append this df (best_n_structures_this_mismatch) to the df best_structure_any_mutation which
        # contains all the best structures per mismatch (regardless of other mutations in the structure) for all genes
        best_structure_any_mutation = best_structure_any_mutation.append(best_n_structures_this_mismatch)

# now that the n best structures for every mismatch have been added to the best_structure_any_mutation df,
# we can sort according to AA index
# To do this, we use the mismatch_of_interest column to create an extra col in best_structure_any_mutation
# containing the AA index of each mismatch (position)
best_structure_any_mutation['aa_index'] = best_structure_any_mutation.mismatch_of_interest.apply(lambda x: int(x[1:-1]))

# now we can sort the df according to the aa_index, mismatch_of_interest and resolution
df_any_mutation = best_structure_any_mutation.sort_values(by=['gene_name', 'aa_index', 'mismatch_of_interest', 'resolution'])
df_any_mutation.reset_index(drop = True, inplace = True)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Add ClinVAr Annotations to dfs
# =======================
# We now have all our dfs ready and can add the respective ClinVar annotations using 
# the functions we defined earlier (to all dfs but the df_wt which won't have any mutations/annotations):
# first we add empty columns to best structure dfs to populate with clinvar data using the function from above:
df_SAV = add_empty_cols(df_SAV)
df_unique_combi = add_empty_cols(df_unique_combi)
df_any_mutation = add_empty_cols(df_any_mutation)

# next, we add clinvar annotations to the new dfs using the function we defined earlier
add_clinvar_annotations(df_SAV)
add_clinvar_annotations(df_unique_combi)
add_clinvar_annotations(df_any_mutation)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Add Validation Report Data from PDB API to dfs
# ===================================
# next, we also want to download additional information (structure validation report) via the PDB API
# we can use the function get_validation_report(pdb_id) which we defined above
# we create an empty df to populate with all the validation report metrics for all the PDB IDs
all_val_reports = pd.DataFrame()

# get a list of all pdb ids in the 4 dfs
# make a set of the list; only contains unique pdb ids
all_pdb_ids = set(list(df_wt.structure_id) + list(df_SAV.structure_id) + list(df_unique_combi.structure_id) + list(df_any_mutation.structure_id))

# now we can iterate over the set with all unique pdb id's to populate the all_val_reports df
for pdb_id in all_pdb_ids:
    validation_report = get_validation_report(pdb_id)
    all_val_reports = all_val_reports.append(validation_report, ignore_index=True)
    
# we drop the column 'resolution' from the all_val_reports, because we don't need this / we already have this info in our dfs (extracted from mmCif files)
all_val_reports.drop(columns='resolution', inplace=True)

# now we can merge the 4 dfs with the all_val_reports df on pdb_id
df_wt = pd.merge(df_wt, all_val_reports, how='left', on='structure_id')
df_SAV = pd.merge(df_SAV, all_val_reports, how='left', on='structure_id')
df_unique_combi = pd.merge(df_unique_combi, all_val_reports, how='left', on='structure_id')
df_any_mutation =  pd.merge(df_any_mutation, all_val_reports, how='left', on='structure_id')

# we have the following dfs with the following cols (see paragraph below for legend)
# df_wt:                              52 cols 
# df_SAV:                           65 cols, incl. (2), (3), (4), (5)
# df_unique_combi:             63 cols, incl (1), (4)
# df_any_mutation:             66 cols, incl (1), (2), (3), (4), (5)

# Legend: some of the dfs have additional columns compared to the df_wt, including:
#     (1) all_mismatches (listing mismatches AND close mismatches):
#     (2) mismatch_of_interest (contains the one mismatch for which we get the best structures): df_any_mutation, df_SAV
#     (3) aa_index (amino acid position of mismatch of interest; integer): df_any_mutation, df_SAV
#     (4) 10 columns for data from ClinVar: 'accession', 'title',
#           'variant_type', 'protein_change', 'aliases', 'clinical significance',
#           'last_evaluated', 'review_status', 'associated_traits',
#           'dbs_and_accessions'
#     (5) mismatch_solved_in_crystal_structure (indicating if mismatch has been solved in structure for each chain)

# For reference: these are all the columns which are potentially in the dfs
# all_cols = [col for col in df_any_mutation.columns]
# print(all_cols)
# ['gene_name', 'structure_id', 'chain_name', 'uniprot_id', 'description', 'species', 'description_ex',
#  'sequence', 'alignment_length', 'hsp_number', 'hsp_length', 'score', 'bit', 'e-value', 'similarity',
#  'mismatches_incl_gaps', 'mismatches_excl_gaps', 'close_mismatches', 'mismatch_substitutions',
#  'close_mismatch_substitutions', 'gaps', 'gaps_in_query', 'gaps_in_sbct', 'gaps_start_pos_query',
#  'gaps_length_query', 'query_del_sbcjt_ins', 'gaps_start_pos_sbjct', 'gaps_length_sbjct',
#  'sbjct_del_query_ins', 'query_sequence', 'sbjct_sequence', 'match_sequence', 'unsolved_residues_in_structure',
#  'n_unsolved_residues', 'percent_unsolved_residues', 'resolution', 'structure_method', 'deposition_date',
#  'structure_name', 'classification', 'all_mismatches', 'mismatch_of_interest', 'mismatch_solved_in_crystal_structure',
#  'aa_index', 'accession', 'title', 'variant_type', 'protein_change', 'aliases', 'clinical significance', 'last_evaluated',
#  'review_status', 'associated_traits', 'dbs_and_accessions', 'r_free', 'r_work', 'r_observed', 'clashscore',
#  'ramachandran_outliers', 'sidechain_outliers', 'RSRZ_outliers', 'polymer_composition', 'oligomeric_state',
#  'symmetry', 'stoichometry', 'full_validation_report']

# first we make a list of columns of interest in the order we want them in the final output files
# this is our list of columns we'd like to have included in the output whenever they are available
output_cols = ['gene_name', 'aa_index', 'mismatch_of_interest', # basic info
               # structure info
               'structure_id', 'structure_name', 'resolution', 'r_free', 'r_work', 'r_observed', 'clashscore',
               'ramachandran_outliers', 'sidechain_outliers', 'RSRZ_outliers', 'polymer_composition', 'oligomeric_state',
               'symmetry', 'stoichometry', 'structure_method', 'mismatch_solved_in_crystal_structure',
               'n_unsolved_residues', 'percent_unsolved_residues', 'unsolved_residues_in_structure',
               'deposition_date', 'classification', 'full_validation_report',
               # blastp info
               'chain_name', 'alignment_length', 'hsp_length', 'score', 'bit', 'e-value', 'similarity', 'mismatch_substitutions',
               'close_mismatch_substitutions', 'gaps',
               # ClinVar info
               'accession', 'title', 'variant_type', 'protein_change', 'aliases', 'clinical significance', 'last_evaluated',
               'review_status', 'associated_traits', 'dbs_and_accessions']

# we need a slightly different order for the unique_combi output!
output_cols_for_unique_combi = ['gene_name', 'mismatch_substitutions', 'close_mismatch_substitutions', # basic info
                # structure info
                'structure_id', 'structure_name', 'resolution', 'r_free', 'r_work', 'r_observed', 'clashscore',
                'ramachandran_outliers', 'sidechain_outliers', 'RSRZ_outliers', 'polymer_composition', 'oligomeric_state',
                'symmetry', 'stoichometry', 'structure_method',  'n_unsolved_residues', 'percent_unsolved_residues',
                'unsolved_residues_in_structure', 'deposition_date', 'classification', 'full_validation_report',
                # blastp info
                'chain_name', 'alignment_length',  'hsp_length', 'score', 'bit', 'e-value', 'similarity', 'gaps',
                # ClinVar info
                'accession', 'title', 'variant_type', 'protein_change', 'aliases', 'clinical significance',
                'last_evaluated', 'review_status', 'associated_traits', 'dbs_and_accessions']

# we use the output_cols list to filter all our dfs and get the columns in the right order before writing to csv
output_wt = df_wt[[col for col in output_cols if col in df_wt.columns]]
output_SAV = df_SAV[[col for col in output_cols if col in df_SAV.columns]]
output_unique_combi = df_unique_combi[[col for col in output_cols_for_unique_combi if col in df_unique_combi.columns]]
output_any_mutation = df_any_mutation[[col for col in output_cols if col in df_any_mutation.columns]]

# Write output files:
# ================
# now we write the different dfs to csv files containing information on all genes to the results folder
print('\nWriting final output for all genes:')
print('    >>> writing csv file called 07_wildtype_structures.csv')
print('            (contains all wildtype structures for all genes')
output_wt.to_csv(f'{results_dir}/07_wildtype_structures.csv', index=False)

print('    >>> writing csv file called 07_best_structures_per_SAV.csv')
print('            (contains the best structures for each SAV for all genes')
output_SAV.to_csv(f'{results_dir}/07_best_structures_per_SAV.csv', index=False)

print('    >>> writing csv file called 07_best_structures_all_unique_combinations.csv')
print('            (contains the best structures for each unique mismatch combination for all genes')
output_unique_combi.to_csv(f'{results_dir}/07_best_structures_all_unique_combinations.csv', index=False)

print('    >>> writing csv file called 07_best_structures_any_mutation.csv')
print('            (contains the best structures for any mutation regardless of other mutations in the same structure for all genes')
output_any_mutation.to_csv(f'{results_dir}/07_best_structures_any_mutation.csv', index=False)

# and finally we write the gene specific outputs in the respective gene folders
# we can loop over all the gene names in the structure_info df (contains all the genes for which structures are available)
for gene in structure_info.gene_name.unique():
    # we get the name of the gene folder of the current gene like so
    try:
        gene_folder = [name for name in os.listdir(results_dir) if name.startswith(f'{gene}_')][0]
    # in case no folder for this gene exists, we skip writing out gene-specific outputs altogether
    except IndexError:
        continue
    
    # get a slice of all the output dfs for just this gene:
    output_slice_wt = output_wt[output_wt.gene_name == gene]
    output_slice_SAV = output_SAV[output_SAV.gene_name == gene]
    output_slice_unique_combi = output_unique_combi[output_unique_combi.gene_name == gene]
    output_slice_any_mutation = output_any_mutation[output_any_mutation.gene_name == gene]    
        
    # and we write the output slice dfs to csv files in the gene folder
    output_slice_wt.to_csv(f'{results_dir}/{gene_folder}/{gene}_07_wildtype_structures.csv', index=False)
    output_slice_SAV.to_csv(f'{results_dir}/{gene_folder}/{gene}_07_best_structures_per_SAV.csv', index=False)
    output_slice_unique_combi.to_csv(f'{results_dir}/{gene_folder}/{gene}_07_best_structures_all_unique_combinations.csv', index=False)
    output_slice_any_mutation.to_csv(f'{results_dir}/{gene_folder}/{gene}_07_best_structures_any_mutation.csv', index=False)


# change back to target directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Extraced best structures for a total of {len(unique_genes)} genes.')
print(f'    hsp coverage threshold: {hsp_coverage}')
print(f'    relative sequence length: {relative_sequence_length}')

print('\nThe following files have been created for each gene and stored in the respective gene folder:')
print('   o      GENENAME_07_best_structures_per_SAV.csv')
print('            (lists best structures for each point mutation (one mutation per structure) in this gene)')    
print('   o      GENENAME_07_best_structures_all_unique_combinations.csv')
print('            (lists best structures for all unique mismatch combinations for this gene)')
print('   o      GENENAME_07_best_structures_any_mutation.csv')
print('            (lists best structures for any mismatch regardless of other mismatches in this gene)')
print('   o      GENENAME_07_wildtype_structures.csv')
print('            (lists all available WT structures for this gene)')

print('\nThe following files have been created and stored in the Results folder:')
print('   o      07_best_structures_per_SAV.csv')
print('            (lists best structures for each point mutation (one mutation per structure) in all genes)')    
print('   o      07_best_structures_all_unique_combinations.csv')
print('            (lists best structures for all unique mismatch combinations for all genes)')
print('   o      07_best_structures_any_mutation.csv')
print('            (lists best structures for any mismatch regardless of other mismatches in all genes)')
print('   o      07_wildtype_structures.csv')
print('            (lists all available WT structures for all genes)\n')

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