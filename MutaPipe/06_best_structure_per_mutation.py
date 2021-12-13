# Script to get best structure per mutation
# This script takes the following csv files as input:
#       - 05_all_info.csv
#       - 02_structure_info.csv
# it will perform the following:
#      - combine the two dfs (according to PDBid)
#      - filter out sequences which are shorter than a given percentage of the reference sequence (set variable relative_sequence_length)
#      - filter out sequences whose best hsp covers less than a given percentage of the reference sequence (set variable hsp_coverage)
#      - sort/filter the df in order to get:
#               - the best structure (best resolution) for all point mutations (structures with only this one mutation and no other mutations)
#               - the best structure (best resolution) for all unique combinations of mutations available in the pdb
#               - the best structure (best resolution) for any specific mutation, regardless of other mutations in the same structure
# and outputs the following files:
#      - In each respective gene folder:
#               - GENENAME_06_best_structure_per_point_mutation.csv
#                 lists best structure for each point mutation (one mutation per structure) in this gene
#               - GENENAME_06_best_structure_all_unique_combinations.csv
#                 lists best structure for all unique mismatch combinations for this gene
#               - GENENAME_06_best_structure_any_mutation.csv
#                 lists best structure for any mismatch regardless of other mismatches in this gene
#      - In the results folder:
#               - 06_best_structure_per_point_mutation.csv
#                  lists best structure for each point mutation (one mutation per structure) in all genes
#               - 06_best_structure_all_unique_combinations.csv
#                  lists best structure for all unique mismatch combinations for all genes
#               - 06_best_structure_any_mutation.csv
#                  lists best structure for any mismatch regardless of other mismatches in all genes

# ===========================================================================================

# Set up
import pandas as pd
import os
import ast
import sys
import argparse
from datetime import datetime

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
                                            
                                            
# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****    Script to get best structure per mutation. This script takes the following csv files as input: 
05_all_info.csv and 02_structure_info.csv. It will perform the following: 
1. combine the two dfs (according to PDBid) 
2. filter out sequences which are shorter than a given percentage of the reference sequence (set variable relative_sequence_length) 
3. filter out sequences whose best hsp covers less than a given percentage of the reference sequence (set variable hsp_coverage) 
4. sort and filter the df in order to get: 
(a) the best structure (best resolution) for all point mutations (structures with only this one mutation and no other mutations) 
(b) the best structure (best resolution) for all unique combinations of mutations available in the pdb 
(c) the best structure (best resolution) for any specific mutation, regardless of other mutations in the same structure 
5. and outputs the following files: 
In each respective gene folder: 
(a) GENENAME_06_best_structure_per_point_mutation.csv (lists best structure for each point mutation [one mutation per structure] in this gene) 
(b) GENENAME_06_best_structure_all_unique_combinations.csv (lists best structure for all unique mismatch combinations for this gene)  
(c) GENENAME_06_best_structure_any_mutation.csv (lists best structure for any mismatch regardless of other mismatches in this gene) 
In the results folder: 
(a) 06_best_structure_per_point_mutation.csv (lists best structure for each point mutation [one mutation per structure] in all genes) 
(b) 06_best_structure_all_unique_combinations.csv (lists best structure for all unique mismatch combinations for all genes) 
(c) 06_best_structure_any_mutation.csv (lists best structure for any mismatch regardless of other mismatches in all genes)    ***""") 

ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')
ap.add_argument("-rsl", "--relative_sequence_length", type=restricted_float, required = False, help=f'filter out sequences shorter than a given percentage of the reference sequence, default = {str(relative_sequence_length)}')
ap.add_argument("-cov", "--hsp_coverage", type=restricted_float, required = False, help=f'filter out sequences whose best hsp covers less than a given percentage of the reference sequence, default = {str(hsp_coverage)}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
relative_sequence_length = relative_sequence_length if args["relative_sequence_length"] == None else args["relative_sequence_length"]
hsp_coverage  = hsp_coverage if args["hsp_coverage"]   == None else args["hsp_coverage"]
# ----------------------------------------------------------------------------------------------------------------------------------

# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# ----------------------------------------------------------------------------------------------------------------------------------

#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_00.txt', 'w') as search_log:
        search_log.write(f'Search log for 00_search_pdb.py\n\n')
    sys.stdout = open(f'{results_dir}/search_log_00.txt', 'a')

# store current date and time in an object and print to console / write to log file
start_time = datetime.now()
print(f'start: {start_time}\n')


# ----------------------------------------------------------------------------------------------------------------------------------

#  read in data
# df = pd.read_csv(f'{results_dir}/07_blast_two_sequences.csv')        #this csv file OR uncomment next line of code
# df = pd.read_csv(f'{results_dir}/06_blast_fasta.csv')
all_info = pd.read_csv(f'{results_dir}/05_all_info.csv')
resolutions = pd.read_csv(f'{results_dir}/02_all_resolutions.csv')
structure_info = pd.read_csv(f'{results_dir}/02_structure_info.csv')

# the blast_two_sequences file does not include the resolutions of the structures, so we read in the resolutions from a previous output
# we combine the two dataframes
structure_info.rename(columns={'gene': 'gene_name'}, inplace=True)
df = all_info.merge(structure_info, how='left', on=['gene_name', 'structure_id'])

# filter out structures which are shorter a given percentage of the reference sequence
# --> relative_sequence_length, e.g. 50% (0.5) or 10% (0.1) of the reference sequence
# ==================================================================================
print(f'Relative sequence length threshold set to {relative_sequence_length}:')
print(f'        only strucures with a sequence length of at least {relative_sequence_length*100}%  of the reference sequence will be included in output\n')
df = df[df.sequence.apply(lambda x: len(x)) > relative_sequence_length]

# filter out structures which cover less than hsp_coverage, e.g. 50% (0.5) or 10% (0.1) of the reference sequence
# ==================================================================================
print(f'Hsp coverage threshold set to {hsp_coverage}:')
print(f'        only strucures with an hsp covering at least {hsp_coverage*100}%  of the reference sequence will be included in output\n')
df = df[df.hsp_length / df.alignment_length > hsp_coverage]


# the table lists the best hsp per structure (in case there are multiple hsps for this strucure's alignment with the reference structure)
# we want to create a filter to exclude very short sequences which only cover a small part of the reference sequence
# the alignment_length in the df is always the same as the length of the reference sequence
# so by dividing the hsp_length / alignment_length, we get the a percentage for the coverage of this hsp for the reference sequence
# (e.g. query covers 75 % of reference sequence)
# with the filter we can exclude all hsps which are below a certain coverage
# we set this threshold to 10% for now
# note: for certain proteins, e.g. FUS, the reference sequence is very long (FUS: 526 AAs) and the hsps of sequences in the structures are much shorter
# for these genes a lower coverage will have to be set if no structures get identified otherwise
# (this applies in all cases where it's difficult to crystallise the entire protein and hence only small fragments are available, I would assume)


# first we create a df that contains the best structure for all point mutations (only one mismatch/close_mismatch)
#  =================================================================================
single_mutation = df[df.close_mismatches + df.mismatches_excl_gaps == 1]

# we order the single mutation df first by gene name, then by normal mismatches, by close mismatches, and finally by resolution
single_mut_sorted = single_mutation.sort_values(by=['gene_name', 'mismatch_substitutions', 'close_mismatch_substitutions', 'resolution'])

# Filter rows
# get rid of all but the best structure per mutation
# Finally we drop all rows but the top/first row (best resolution) for each gene and mutation
best_structure_single_mut = single_mut_sorted.drop_duplicates(subset=['gene_name', 'mismatch_substitutions', 'close_mismatch_substitutions'], keep="first")


# in order to get the aa index for each mutation so we are able to sort the df accordingly, we do the following:
# make mismatch of interest col in df best_structure_single_mut
# first we make a column called mismiatch_of_interest, so we can do the same as above afterwards
best_structure_single_mut.loc[:,'mismatch_of_interest'] = None
# to fill the column with values, we loop over the df (couldn't figure out how else to do it)
for index, row in best_structure_single_mut.iterrows():
    close_mis = ast.literal_eval(row.close_mismatch_substitutions)
    mis = ast.literal_eval(row.mismatch_substitutions)
    all_mismatches = mis +close_mis
    best_structure_single_mut.loc[index, 'mismatch_of_interest'] = all_mismatches[0]
    
# now we use mismatch_of_interest column to create extra col in best_structure_single_mut
# containing the AA index of each mismatch (position)
best_structure_single_mut.loc[:,'aa_index'] = best_structure_single_mut.mismatch_of_interest.apply(lambda x: int(x[1:-1]))

# now we can sort the df according to the aa_index
best_structure_single_mut = best_structure_single_mut.sort_values(by=['gene_name', 'aa_index'])

# write only columns of interest to ouput for now:
output = best_structure_single_mut[['gene_name', 'aa_index', 'mismatch_of_interest', 'structure_id', 'structure_name', 'chain_name', 'resolution', 'mismatch_substitutions', 'close_mismatch_substitutions', 'unsolved_residues_in_structure', 'structure_method', 'deposition_date', 'classification']]

# write final df to file
print('>>> writing csv file called 06_best_structure_per_point_mutation.csv containing the best structure for each point mutation for all genes\n')
output.to_csv(f'{results_dir}/06_best_structure_per_point_mutation.csv', index=False)


# get the best structure for any unique combination of mutations per gene
# =====================================================

# in order to do that we first create a new column in the df which contains a list of all the mismatch substitutions
# (close mismatches and normal mismatches)
# we make two lists of the mismatch and close_mismatch column and add them to each other
# we also check if the variable is a string. missing values are floats, so they get ignored as we write an empty list instead
df['all_mismatches'] = df.mismatch_substitutions.apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else []) + df.close_mismatch_substitutions.apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else [])

# create a list of all unique gene names, so we can loop over a slice of the df containing only structures of this gene
# we do this because the mismatches for one gene are comparable, but not for multiple genes (e.g. A5E is not the same mutation if it is in two differenct proteins and cannot be summed up in one row)
unique_genes = df.gene_name.unique()

# create empty df to populate with best structure per unique combination:
best_structure_unique_combis = pd.DataFrame(columns=df.columns)

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
    # we can then filter again in order to get the structure with the highest resolution for each combination!

    # we do this by looping over the unique values, extracting the matching rows and then selecting the best structure
    # (sorting by structure_id, resolution) and then dropping all duplicates but keeping the first row (=best resolution)
    
    # we also create an empty df to populate with the best structures per unique combination for this gene only
    # we will write this output to the respective gene folder
    combis_per_gene = pd.DataFrame(columns = df.columns)
    # initiate for loop
#     print(f'>>> getting best structure for each unique mismatch combination for {gene}')

    for combination in unique_mismatch_combinations:
        # we extract a slice of the df only containing the rows with this mismatch combination
        df_this_combi = df_gene[df_gene.all_mismatches == combination]
        
        # sort by resolution
        df_this_combi_sorted = df_this_combi.sort_values(by='resolution')
        
        # drop all but the best resolution (keep first row)
        best_structure_this_combi = df_this_combi_sorted.drop_duplicates(subset='all_mismatches', keep="first")
        
        # we now append the best_structure_this_combi to the best_structure_unique_combis df
        best_structure_unique_combis = best_structure_unique_combis.append(best_structure_this_combi)
        
        # we also append it to the gene-specific df combis_per_gene
        combis_per_gene = combis_per_gene.append(best_structure_this_combi)
        
    # now that we have the best strucutre for all unique mismatch combinations for this gene,
    # we want to write the gene-specific df to a csv file in the respective gene folder.
    # first we need to identify the right gene folder.
    # it's in the results folder and starts with the string '{gene}_'
    gene_folder = [name for name in os.listdir(results_dir) if name.startswith(f'{gene}_')][0]
    
    # extrac the relevant columns and write to csv file
    print('    >>> writing csv file called {gene}_06_best_structure_per_unique_combination.csv')
    print(f'            (contains the best structure for each unique mismatch combination for {gene}')
    combis_per_gene_output = combis_per_gene[['gene_name', 'mismatch_substitutions', 'close_mismatch_substitutions', 'structure_id', 'structure_name', 'chain_name', 'resolution', 'unsolved_residues_in_structure', 'structure_method', 'deposition_date', 'classification']]
    combis_per_gene_output.to_csv(f'{results_dir}/{gene_folder}/{gene}_06_best_structure_per_unique_combination.csv', index = False)
    
    # we also want to store a csv file containing only best structure per the point mutation for this gene:
    # we previously stored the df containing only the relevant columns for ALL genes in a df called 'output' (see above)
    # we get a slice of the output df for this particular gene and then write it to a csv file in the respective folder
    # the correct folder is results_dir/gene_folder
    print('    >>> writing csv file called {gene}_06_best_structure_per_unique_combination.csv')
    print(f'            (contains the best structure for each point mutation for {gene}')
    output_gene = output[output.gene_name == gene]
    output_gene.to_csv(f'{results_dir}/{gene_folder}/{gene}_06_best_structure_per_point_mutation.csv', index=False)
    

    # get the best structure for any mutation regardless of whether there are other mutations (also per gene, still in the loop)
    # =========================================================================================

    # reconvert the all_mismatches column in df_gene from string to list
    df_gene['all_mismatches'] = df_gene.all_mismatches.apply(lambda x: ast.literal_eval(x))
    
    # we get a list of all unique mismatches (not mismatch combinations, but single point mutations) for the current gene
    unique_mismatches = []
    
    # we loop over the all_mismatches columns,
    # each row contains a list with all mismatches in this structure
    for mismatches in df_gene.all_mismatches:
        # we loop over all the mismatches in this structure
        for mismatch in mismatches:
            # we append this mismatch to unique_mismatches if it's not already there
            if mismatch not in unique_mismatches:
                unique_mismatches.append(mismatch)
    
    # we create an empty df to populate with the best structure for each unique mismatch for the current gene
    mutations_per_gene = pd.DataFrame(columns = df.columns)
    # we add a empty column where we can store each mismatch so it's clear for which mismatch this is the best strucure
    # (as there might be multiple mismatches in the structure overall)
    mutations_per_gene['mismatch_of_interest'] = None

    # we have to convert the all_mismatches column back to a string (it's currently a list)
    df_gene['all_mismatches'] = df_gene.all_mismatches.apply(lambda x: str(x))

    # now we loop over the list of unique mismatches for this gene and get the best structure for each of them
#     print(f'>>> getting best structure for any mutation regardless of other mutations for {gene}')
    for mismatch in unique_mismatches:
        # first we get a df of all structures which contain this mismatch (regardless of other mismatches in the structure)
        this_mismatch = df_gene[df_gene.all_mismatches.str.contains(mismatch)]
        # now we sort the this_mismatch df to get only the best structure (highest resolution)
        # sort by  resolution only
        this_mismatch_sorted = this_mismatch.sort_values(by=['resolution'])
        
        # drop all but the best resolution (keep first row)
        # this df is only for one mismatch, so it only has one row after we drop all but the best structure!
        best_structure_this_mismatch = this_mismatch_sorted.drop_duplicates(subset=['gene_name', ], keep="first").copy()
        
        # we add information on this mismatch to the new column we create earlier
        best_structure_this_mismatch['mismatch_of_interest'] = mismatch
        
        # we append this df (best_structure_this_mismatch) to the df mutations_per_gene which
        # contains all the best structures per mismatch (regardless of other mutations in the structure)
        mutations_per_gene = mutations_per_gene.append(best_structure_this_mismatch)
        
    # now that we looped over all the unique mismatches, we extrac the relevant columns from the df and write to csv file for this gene
    print('    >>> writing csv file called {gene}_06_best_structure_any_mutation.csv')
    print(f'            (contains the best structure for any mutation regardless of other mutations for {gene}\n')
    mutations_per_gene_output = mutations_per_gene[['gene_name', 'mismatch_of_interest', 'structure_id', 'structure_name', 'chain_name', 'resolution', 'mismatch_substitutions', 'close_mismatch_substitutions', 'unsolved_residues_in_structure', 'structure_method', 'deposition_date', 'classification']]
    mutations_per_gene_output.to_csv(f'{results_dir}/{gene_folder}/{gene}_06_best_structure_any_mutation.csv', index = False)
    
    # we also append the mutations_per_gene df to the best_structure_any_mutation df
    best_structure_any_mutation = best_structure_any_mutation.append(mutations_per_gene)

#  sort according to AA index
# use mismatch_of_interest column to create extra col in best_structure_any_mutation
# containing the AA index of each mismatch (position)
best_structure_any_mutation['aa_index'] = best_structure_any_mutation.mismatch_of_interest.apply(lambda x: int(x[1:-1]))

# now we can sort the df according to the aa_index
best_structure_any_mutation = best_structure_any_mutation.sort_values(by=['gene_name', 'aa_index'])


# now we write the csv files containing information on all genes to the results folder

# extract relevant columns for output file
output1 = best_structure_unique_combis[['gene_name', 'mismatch_substitutions', 'close_mismatch_substitutions', 'structure_id', 'structure_name', 'chain_name', 'resolution', 'unsolved_residues_in_structure', 'structure_method', 'deposition_date', 'classification']]
output2 = best_structure_any_mutation[['gene_name', 'aa_index', 'mismatch_of_interest', 'structure_id', 'structure_name', 'chain_name', 'resolution', 'mismatch_substitutions', 'close_mismatch_substitutions', 'unsolved_residues_in_structure', 'structure_method', 'deposition_date', 'classification']]
# write final dfs to file
print('Writing final output for all genes')
print('    >>> writing csv file called 06_best_structure_all_unique_combinations.csv')
print(f'            (contains the best structure for each unique mismatch combination for all genes')
output1.to_csv(f'{results_dir}/06_best_structure_all_unique_combinations.csv', index=False)
print('    >>> writing csv file called 06_best_structure_all_unique_combinations.csv')
print(f'            (contains the best structure for any mutation regardless of other mutations for all genes')
output2.to_csv(f'{results_dir}/06_best_structure_any_mutation.csv', index=False)    


print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Extraced best structures for a total of {len(unique_genes)} genes.')
print(f'    hsp coverage threshold: {hsp_coverage}')
print(f'    relative sequence length: {relative_sequence_length}')

print('\nThe following files have been created for each gene and stored in the respective gene folder:')
print('   o      GENENAME_06_best_structure_per_point_mutation.csv')
print('            (lists best structure for each point mutation (one mutation per structure) in this gene)')    
print('   o      GENENAME_06_best_structure_all_unique_combinations.csv')
print('            (lists best structure for all unique mismatch combinations for this gene)')
print('   o      GENENAME_06_best_structure_any_mutation.csv')
print('            (lists best structure for any mismatch regardless of other mismatches in this gene)')    

print('\nThe following files have been created and stored in the Results folder:')
print('   o      06_best_structure_per_point_mutation.csv')
print('            (lists best structure for each point mutation (one mutation per structure) in all genes)')    
print('   o      06_best_structure_all_unique_combinations.csv')
print('            (lists best structure for all unique mismatch combinations for all genes)')
print('   o      06_best_structure_any_mutation.csv')
print('            (lists best structure for any mismatch regardless of other mismatches in all genes)\n')


# store current date and time in an object and print to console / write to log file
end_time = datetime.now()
print(f'start: {start_time}\n')
print(f'end: {end_time}\n\n')

# close search log
if create_search_log == True:
    sys.stdout.close()


