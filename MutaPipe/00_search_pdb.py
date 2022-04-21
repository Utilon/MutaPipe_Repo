# This script takes a list of genes in txt. format as input and performs the following:
#       - searches the pdb for all structures associated with each gene name (in Homo Sapiens)
#       - outputs a csv file called 00_search_overview_PDBids.csv containing all gene names and corresponding PDB structure id's if available
#       - outputs a csv file called 00_search_overview_availability.csv containing a boolean value for each gene to indicate whether there are any structures available or not


# ----------------------------------------------------------------------------------------------------------------------------------
# Set up workspace
import os
import pandas as pd
import sys
import argparse
import requests
from datetime import datetime

# get this script's name:
script_name = os.path.basename(__file__)
# ----------------------------------------------------------------------------------------------------------------------------------

#  read in gene_list from textfile
with open('genes.txt', 'r') as f:
    gene_list = f.read().split(' ')
    # remove empty entries
    if '' in gene_list:
        gene_list.remove('')

# set default values for arguments we want to implement
# we have to do this here if we want to print the default values in the help message
# specify other search terms and search operator:
species_name = "Homo sapiens" # script will search for exact match 
log_operator = "and"
all_hits = True          # if set to "false": output limited to 10 pdb ids per gene
target_directory = os.getcwd()   # set directory to create Results folder with new folders with pdb files for each gene
                                                # Default: set to current working directory (where this script is saved)
create_search_log = False      # will create a file called search_log.txt with console output if set to True,
                                            # prints to console if set to False.
                                            
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

# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****   This script takes a list of genes in txt. format as input and performs the following:
1. searches the pdb for all structures associated with each gene name (in Homo Sapiens)
2. outputs a csv file called 00_search_overview_PDBids.csv containing all gene names and corresponding PDB structure id's if available
3. outputs a csv file called 00_search_overview_availability.csv containing a boolean value for each gene to indicate whether there are any structures available or not   ***""")

ap.add_argument('-g','--genes', nargs='+', required=False, help=f'Specify genes for which to search pdb structures, default = {gene_list}; to pass a file containing all genes use -g $(cat filename)')
ap.add_argument("-o", "--organism", required = False, help=f'Specify species for which to search pdb structures, default = {species_name}')
ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'Write output to .log file in current directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'Specify target directory, default = {target_directory}')
ap.add_argument("-a", "--all", type=str2bool, required = False, help=f'Retrieve all (True) vs max. 10 pdb IDs per gene (False), default = {str(all_hits)}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
gene_list = gene_list if args["genes"] == None else args["genes"]
species_name = species_name if args["organism"] == None else args["organism"]
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
all_hits = all_hits if args["all"] == None else args["all"]

# ----------------------------------------------------------------------------------------------------------------------------------

# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# create Results folder if it doesn't already exist
if not os.path.exists(results_dir):
        os.makedirs(results_dir)

# ----------------------------------------------------------------------------------------------------------------------------------

#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_00.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name} with genes {gene_list}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_00.txt', 'a')

# print nice title
print('===============================================================================')
print('*****    Searching the Protein Data Bank for Structures Associated with Input Genes    *****')
print('===============================================================================\n')

# print script name to console/log file
print(f'script name: {script_name}')

# store current date and time in an object and print to console / write to log file
start_time = datetime.now()
print(f'start: {start_time}\n')

# ----------------------------------------------------------------------------------------------------------------------------------

# set up search request
base_url = "https://search.rcsb.org/rcsbsearch/v1/query"

# define search function to get data from API
def search_pdb(query_url):
    """
    This function will make a GET call to the PDB API
    using the defined query url
    
    :param query_url: string
    :return: Dict or None
    """
    
    # Make a GET call to the API URL
    get_request = requests.get(url=query_url)
    
    if get_request.status_code == 200:
        # If there is data returned (with HTML status code 200)
        # then return the data in JSON format
        return get_request.json()
    else:
        # If there is no data, print status code and response, return None
        print(f'            No data retrieved for {gene_name}; status code: {get_request.status_code} {get_request.text}\n')
        return None

# ----------------------------------------------------------------------------------------------------------------------------------

# change to results directory:
os.chdir(results_dir)

# Initiate loop over all genes in gene_list to find all PDB id's of available structures
gene_counter = 1
n_genes = len(gene_list)
n_pdbs_all_genes = 0
genes_data_available = []
genes_no_data_available = []
search_overview = pd.DataFrame(columns=['gene_name', 'n_available_structures', 'available_structures'])

print(f'\n======================== Searching PDB IDs for {n_genes} inputted genes ========================\n')

for gene in gene_list:
    gene_name = gene
    
    # Required arguments for all nodes in search sent to pdb 
    query_type1 = "group"     # "terminal" or " group"
    
    query_type2 = "terminal"
    query_service2 = "text"     # also possible to search for 'sequence', 'seqmotif', 'structure', 'strucmotif' etc. (see documentation)
    
    query_type3 = "terminal"
    query_service3 = "text"
    
    return_type = "entry"       # 'assembly', 'polymer_entitiy' etc. (see documentation)
    
    # Optional arguments for terminal nodes
    parameter_attribute2 = "rcsb_entity_source_organism.rcsb_gene_name.value" # we want to search for a specific gene
    parameter_operator2 = "exact_match"
    parameter_value2 =  gene_name.replace(" ", "%20") #replace spaces in search term with '%20' to make it url-conform
    
    parameter_attribute3 = "rcsb_entity_source_organism.taxonomy_lineage.name" # we want to search for a species (i.e. Homo Sapiens')
    parameter_operator3 = "exact_match"
    parameter_value3 =  species_name.replace(" ", "%20") #replace spaces in search term with '%20' to make it url-conform
    
    # ----------------------------------------------------------------------------------------------------------------------------------
    
    # Build query string (no packages required)
    # Make json string url-compatible
    
    # Currently this script only works with this exact query format
    # in json format, this would like like this (without the three quotation marks on each side)
    # query_string_json = """{"query":{"type":"group","logical_operator":"and","nodes":[{"type":"terminal","service":"text","parameters":{"operator":"exact_match","value":"Homo sapiens","attribute":"rcsb_entity_source_organism.taxonomy_lineage.name"}},{"type":"terminal","service":"text","parameters":{"operator":"exact_match","value":"TARDBP","attribute":"rcsb_entity_source_organism.rcsb_gene_name.value"}}]},"return_type":"entry"}"""
    # first i have replaced all curly brackets with appropriate url-code '%7B' and '%7D', like so:
    # query_string_replace1 = query_string_json.replace("{", "%7B")
    # query_string_replace2 = query_string_replace1.replace("}", "%7D")
    # then I have manually inserted the variables into query_string_replace2 (that's why only certain parameters are adjustable in this script atm)
    #  query_string_with_vars = f"""%7B"query":%7B"type":"{query_type1}","logical_operator":"{log_operator}","nodes":[%7B"type":"{query_type2}","service":"{query_service2}","parameters":%7B"operator":"{parameter_operator2}","value":"{parameter_value2}","attribute":"{parameter_attribute2}"%7D%7D,%7B"type":"{query_type3}","service":"{query_service3}","parameters":%7B"operator":"{parameter_operator3}","value":"{parameter_value3}","attribute":"{parameter_attribute3}"%7D%7D]%7D,"request_options":%7B"return_all_hits":{all_hits}%7D,"return_type":"{return_type}"%7D"""
    query_string_with_vars = f"""%7B"query":%7B"type":"{query_type1}","logical_operator":"{log_operator}","nodes":[%7B"type":"{query_type2}","service":"{query_service2}","parameters":%7B"operator":"{parameter_operator2}","value":"{parameter_value2}","attribute":"{parameter_attribute2}"%7D%7D,%7B"type":"{query_type3}","service":"{query_service3}","parameters":%7B"operator":"{parameter_operator3}","value":"{parameter_value3}","attribute":"{parameter_attribute3}"%7D%7D]%7D,"request_options":%7B"return_all_hits":{str(all_hits).lower()}%7D,"return_type":"{return_type}"%7D"""
    # now we can go on and replace " and square brackets from the query string with url-code
    query_string_replace3 = query_string_with_vars.replace('"', '%22')
    query_string_replace4 = query_string_replace3.replace("[", "%5B")
    query_string = query_string_replace4.replace("]", "%5D") # we call our final query link 'query_string'
    
    # ----------------------------------------------------------------------------------------------------------------------------------
    
    # Set up search
    query_url = base_url + "?json=" + query_string
    # https://search.rcsb.org/rcsbsearch/v1/query?json={search-request}
    # the curly brackets around the search-request are defined in the query_string
    
    # ----------------------------------------------------------------------------------------------------------------------------------
    
    # perfom search
       
    print(f">>> Searching the PDB for gene {gene_counter} of {n_genes}:             {gene_name}")
    print(f'            Search terms:                            "{gene_name}" {log_operator.upper()} "{species_name}"')
    
    pdb_data = search_pdb(query_url)
    # this is actual data coming from the PDB API
    # as this contains information we don't need (metadata etc.), we isolate the results
    # and then create a list containing all pdb ids ("identifiers")
    try:
        results = pdb_data['result_set']
        genes_data_available.append(gene_name)
    # if we can't retrieve data and get an error instead, we skip this loop and go to the next gene in the list
    except:
        genes_no_data_available.append(gene_name)
        gene_counter += 1
        continue
    
    found_pdbs = []
    for entry in results:
        pdb_id = entry["identifier"]
        found_pdbs.append(pdb_id.lower())
    
    # sort found_pdbs
    found_pdbs.sort()
    print(f'            Number of PDB IDs retrieved :            {len(found_pdbs)} \n')
    
    # add number of found pdb ID's to counter n_pdbs_all_genes
    n_pdbs_all_genes += len(found_pdbs)
    
    # append information to search_overview df to later write to file
    # in order to add a new row to the bottom of the df, we use .loc and specify the index as the current df lenght!
    search_overview.loc[len(search_overview)] = [gene, len(found_pdbs), found_pdbs]
    gene_counter += 1


# write search_overview to csv file
search_overview.to_csv(f'{results_dir}/00_search_overview_PDBids.csv')

# write a csv file that contains a boolean value for each gene to indicate whether there are any structures available or not

# create an empty dataframe
df = pd.DataFrame(columns=['gene_name', 'data_available'])

# populate df with information on availability
for gene in gene_list:
    if gene in genes_data_available:
        df.loc[len(df)] = [gene, True]
    elif gene in genes_no_data_available:
        df.loc[len(df)] = [gene, False]
    else:
        df.loc[len(df)] = [gene, np.nan]
        
df.to_csv(f'{results_dir}/00_search_overview_availability.csv', index = False)

# change back to target directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'    o      A total of {n_pdbs_all_genes} PDB IDs have been found for the inputted gene list.')
print(f'    o      PDB IDs have been found for {len(genes_data_available)} out of {n_genes} gene(s) on your list.')
print(f'    o      No PDB IDs have been found for {len(genes_no_data_available)} gene(s).\n')

print('The following files have been created:')
print('   o      00_search_overview_availability.csv      (contains information on which genes have available structures)')
print('   o      00_search_overview_PDBids.csv            (contains information on which PDB IDs are available per structures)\n\n')



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