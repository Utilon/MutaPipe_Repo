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
(6) output the following files:
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

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
relative_sequence_length = relative_sequence_length if args["relative_sequence_length"] == None else args["relative_sequence_length"]
hsp_coverage  = hsp_coverage if args["hsp_coverage"]   == None else args["hsp_coverage"]
n_best_structures = n_best_structures if args["n_best_structures"] == None else args["n_best_structures"]

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

for index, row in blastp_results.iterrows():
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
    # and attach this information to the blastp_results df
    # first we create an empty dictionary to populate with all chains (keys) and unsolved residues (values) in this structure
    # which have the same sequence and are thus in one and the same row in the blastp_results df
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

    # now that we've looped over all chains which have the same sequence in this structure, we add the unsolved_dict to the blastp_results df
    blastp_results.loc[index, 'unsolved_residues_in_structure'] = str(unsolved_dict)

# now the blastp_results df contains all corresponding unsolved residues!
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

# Filter rows
# in order to get the aa index for each mutation so we are able to sort the df accordingly, we do the following:
# make mismatch of interest col in df single_mut_sorted 
# first we make a column called mismiatch_of_interest
# single_mut_sorted.loc[:,'mismatch_of_interest'] = None # this line does not work with KL 24.6.2022
# try this instead
single_mut_sorted['mismatch_of_interest'] = None  # this works with KL 24.6.2022
# to fill the column with values, we loop over the df (couldn't figure out how else to do it)
for index, row in single_mut_sorted.iterrows():
    close_mis = ast.literal_eval(row.close_mismatch_substitutions)
    mis = ast.literal_eval(row.mismatch_substitutions)
    all_mismatches = mis +close_mis
    single_mut_sorted.loc[index, 'mismatch_of_interest'] = all_mismatches[0]
    
# now we use mismatch_of_interest column to create extra col in single_mut_sorted
# containing the AA index of each mismatch (position)
single_mut_sorted.loc[:,'aa_index'] = single_mut_sorted.mismatch_of_interest.apply(lambda x: int(x[1:-1]))

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
    # if this throws an Error we don't sort the df, and use the other df with all relevant column names
    df_SAV = single_mut_sorted
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
# this df also needs an extra column called mismatch_of_interes
best_structure_any_mutation['mismatch_of_interest'] = None

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
        
        # FILTER n_best_structures        
        # we want to get the first n rows of the this_mismatch_sorted df
        best_n_structures_this_mismatch = this_mismatch_sorted[:n_best_structures]        
        
        # we add information on this mismatch to the new column we created earlier
        best_n_structures_this_mismatch['mismatch_of_interest'] = mismatch
        
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

# we have the following dfs with the following cols (see paragraph below for legend)
# df_wt:                              38 cols
# df_SAV:                           50 cols, incl. (2), (3), (4)
# df_unique_combi:             49 cols, incl (1), (4)
# df_any_mutation:             51 cols, incl (1), (2), (3), (4)

# Legend: some of the dfs have additional columns compared to the df_wt, including:
#     (1) all_mismatches (listing mismatches AND close mismatches):
#     (2) mismatch_of_interest (contains the one mismatch for which we get the best structures): df_any_mutation, df_SAV
#     (3) aa_index (amino acid position of mismatch of interest; integer): df_any_mutation, df_SAV
#     (4) 10 columns for data from ClinVar: 'accession', 'title',
#           'variant_type', 'protein_change', 'aliases', 'clinical significance',
#           'last_evaluated', 'review_status', 'associated_traits',
#           'dbs_and_accessions'

# For reference: these are all the columns which are potentially in the dfs
# all_cols = [col for col in df_unique_combi.columns]
# print(all_cols)
# ['gene_name', 'structure_id', 'chain_name', 'uniprot_id', 'description', 'species',
#  'description_ex', 'sequence', 'alignment_length', 'hsp_number', 'hsp_length', 'score', 'bit', 'e-value', 'similarity',
#  'mismatches_incl_gaps', 'mismatches_excl_gaps', 'close_mismatches', 'mismatch_substitutions',
#  'close_mismatch_substitutions', 'gaps', 'gaps_in_query', 'gaps_in_sbct', 'gaps_start_pos_query',
#  'gaps_length_query', 'query_del_sbcjt_ins', 'gaps_start_pos_sbjct', 'gaps_length_sbjct', 'sbjct_del_query_ins',
#  'query_sequence', 'sbjct_sequence', 'match_sequence', 'unsolved_residues_in_structure', 'resolution',
#  'structure_method', 'deposition_date', 'structure_name', 'classification', 'all_mismatches',
#  'accession', 'title', 'variant_type', 'protein_change', 'aliases', 'clinical significance',
#  'last_evaluated', 'review_status', 'associated_traits', 'dbs_and_accessions']

# first we make a list of columns of interest in the order we want them in the final output files
# this is our list of columns we'd like to have included in the output whenever they are available
output_cols = ['gene_name', 'aa_index', 'mismatch_of_interest', 'structure_id', 'structure_name', 'chain_name', 'resolution',
    'alignment_length',  'hsp_length', 'score', 'bit', 'e-value', 'similarity',
    'mismatch_substitutions', 'close_mismatch_substitutions', 'gaps',
    'structure_method', 'deposition_date', 'classification', 'unsolved_residues_in_structure',
    'accession', 'title', 'variant_type', 'protein_change', 'aliases', 'clinical significance',
    'last_evaluated', 'review_status', 'associated_traits', 'dbs_and_accessions']

# we need a slightly different order for the unique_combi output!
output_cols_for_unique_combi = ['gene_name', 'mismatch_substitutions', 'close_mismatch_substitutions',
    'structure_id', 'structure_name', 'chain_name', 'resolution', 
    'alignment_length',  'hsp_length', 'score', 'bit', 'e-value', 'similarity', 'gaps',
    'structure_method', 'deposition_date', 'classification', 'unsolved_residues_in_structure',
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