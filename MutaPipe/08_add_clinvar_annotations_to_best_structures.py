# Combine ClinVar annotations with best_structure outputs
# ************************************************************************************
# This script takes the following csv files as input
#      - 06_best_structure_all_unique_combinations.csv
#      - 06_best_structure_any_mutation.csv
#      - 06_best_structure_per_point_mutation.csv
#      - 07_b_ClinVar_Annotations.csv
#  and will:
#      - add availalable ClinVar annotations to all three best_structure tables
#      - output the following files:
#          - 08_best_structure_all_unique_combinations.csv (incl. ClinVar annotations)
#          - 08_best_structure_any_mutation.csv (incl. ClinVar annotations)
#          - 08_best_structure_per_point_mutation.csv (incl. ClinVar annotations)
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
  
# set default values for arguments we want to implement
# we have to do this here if we want to print the default values in the help message

create_search_log = False     # will create a file called search_log.txt with console output if set to True,
                                            # prints to console if set to False.

target_directory = os.getcwd()    # set target directory (where Results folder is located)
                                            
                                            
# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****    This script takes the following csv files as input: 
(a) 06_best_structure_all_unique_combinations.csv 
(b) 06_best_structure_any_mutation.csv 
(c) 06_best_structure_per_point_mutation.csv 
(d) 07_b_ClinVar_Annotations.csv 
and will: 
1. add availalable ClinVar annotations to all three best_structure tables 
2. output the following files: 
(a) 08_best_structure_all_unique_combinations.csv (incl. ClinVar annotations) 
(b) 08_best_structure_any_mutation.csv (incl. ClinVar annotations) 
(c) 08_best_structure_per_point_mutation.csv (incl. ClinVar annotations)    ***""")

ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]

# ----------------------------------------------------------------------------------------------------------------------------------

# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# ----------------------------------------------------------------------------------------------------------------------------------

#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_08.txt', 'w') as search_log:
        search_log.write(f'Search log for 08_add_clinvar_annotations_to_best_structures.py\n\n')
    sys.stdout = open(f'{results_dir}/search_log_08.txt', 'a')

# print nice title
print('===============================================================================')
print('*****    Adding Relevant ClinVar Annotations to Best Structure Data    *****')
print('===============================================================================\n')

# print script name to console/log file
print(f'script name: {os.path.basename(__file__)}')


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
    for index, row in df.iterrows():
        gene = row.gene_name
        # get a slice of the clinvar_annotations df with variants for only the current gene
        clinvar_annotations_this_gene = clinvar_annotations[clinvar_annotations.gene == gene]
        
        # we add a try and except statement to use the mismatch_of_interest colum if available
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

# read in data
unique_combi = pd.read_csv(f'{results_dir}/06_best_structure_all_unique_combinations.csv')
any_mutation = pd.read_csv(f'{results_dir}/06_best_structure_any_mutation.csv')
point_mutations = pd.read_csv(f'{results_dir}/06_best_structure_per_point_mutation.csv')
clinvar_annotations = pd.read_csv(f'{results_dir}/07_b_ClinVar_Annotations.csv')

# first we add empty columns to best structure dfs to populate with clinvar data using the function from above:
new_unique_combi = add_empty_cols(unique_combi)
new_any_mutation = add_empty_cols(any_mutation)
new_point_mutations = add_empty_cols(point_mutations)

# next, we add clinvar annotations to the new dfs using the function we defined earlier
add_clinvar_annotations(new_point_mutations)
add_clinvar_annotations(new_any_mutation)
add_clinvar_annotations(new_unique_combi)

# now we can write the new dfs to csv files
new_unique_combi.to_csv(f'{results_dir}/08_best_structure_all_unique_combinations.csv', index = False)
new_any_mutation.to_csv(f'{results_dir}/08_best_structure_any_mutation.csv', index = False)
new_point_mutations.to_csv(f'{results_dir}/08_best_structure_per_point_mutation.csv', index = False)

print('\n============================== Summary ================================================\n')
print('Complete! \n    Added all available ClinVar annotations.\n')

print('The following files have been created and stored in the Results folder:')
print('   o      08_best_structure_all_unique_combinations.csv')
print('            (lists best structure for all unique mismatch combinations for all genes, incl. ClinVar annotations)')
print('   o      08_best_structure_any_mutation.csv')
print('            (lists best structure for any mismatch regardless of other mismatches in all genes, incl. ClinVar annotations)')
print('   o      08_best_structure_per_point_mutation.txt')
print('            (lists best structure for each point mutation (one mutation per structure) in all genes, incl. ClinVar annotations)\n')


# print script name to console/log file
print(f'end of script {os.path.basename(__file__)}')

# store current date and time in an object and print to console / write to log file
end_time = datetime.now()
print(f'start: {start_time}')
print(f'end: {end_time}\n\n')
print('........................................................................................................................................................\n\n\n')

# close search log
if create_search_log == True:
    sys.stdout.close()