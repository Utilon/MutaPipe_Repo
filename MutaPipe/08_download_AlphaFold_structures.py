# Script to download AlphaFold predicted structures for MutaPipe input genes
# ==========================================================================================
# Script will use input gene names and get the corresponding UniProt ID from the UniProt API
# The UniProt ID is then used to construct the url to download the AlphaFold predicted structures from the AlphaFold database
# ==========================================================================================
# Set up
# ******
import requests
import pandas as pd
import argparse
import os
from datetime import datetime
import sys

# get this script's name:
script_name = os.path.basename(__file__)
# ----------------------------------------------------------------------------------------------------------------------------------
# Read in gene file for DEFAULT settings
#  read in gene_list from textfile
try:
    with open('genes.txt', 'r') as f:
        genes = f.read().split(' ')
    # remove empty entries
    if '' in genes:
        genes.remove('')
except FileNotFoundError:
    genes = ['DCTN1', 'ERBB4', 'SOD1']
# set default values for arguments we want to implement
# we have to do this here if we want to print the default values in the help message
# specify other search terms and search operator:
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
ap = argparse.ArgumentParser(description="""****   This script takes a list of genes in txt. format as input and performs the following:\n
1. Gets the corresponding UniProt IDs for each gene name (in Homo Sapiens) via the UniProt API\n
2. Downloads all AlphaFold2 predicted structures for the identified UniProt IDs\n
3. Outputs a csv file called 08_AlphaFold_structures indicating download status for each structure   ***""")

ap.add_argument('-g','--genes', nargs='+', required=False, help=f'Specify genes for which to download AlphaFold predicted structures, default = {genes}; to pass a file containing all genes use -g $(cat filename)')
ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'Write output to .log file in current directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'Specify target directory, default = {target_directory}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
genes = genes if args["genes"] == None else args["genes"]
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]

# ----------------------------------------------------------------------------------------------------------------------------------
# We want to write all our Output into the Results directory
results_dir = f'{target_directory}/Results' #define path to results directory

# create Results folder if it doesn't already exist
if not os.path.exists(results_dir):
        os.makedirs(results_dir)
# ----------------------------------------------------------------------------------------------------------------------------------
#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_08.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name} with genes {genes}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_08.txt', 'a')

# print nice title
print('===============================================================================')
print('*****    Downloading AlphaFold2 Predicted Structures for all Input Genes    *****')
print('===============================================================================\n')

# print script name to console/log file
print(f'script name: {script_name}')

# store current date and time in an object and print to console / write to log file
start_time = datetime.now()
print(f'start: {start_time}\n')

print('Input genes: ', genes, '\n')
# ----------------------------------------------------------------------------------------------------------------------------------

# ==========================================================================================
# define functions with GET requests
# **********************************

# Define a function which will retrieve data from UniProt for all input gene names
# see the following links for more info:
# this is basically an example of a search result which we might want to download:
# (click on download to choose format + settings and click generate url for API)
# https://www.uniprot.org/uniprotkb?query=organism_name%3A%22homo%20sapiens%22%20AND%20%28gene_exact%3Abraf%20OR%20gene_exact%3Abrca1%20OR%20gene_exact%3Abrca2%20OR%20gene_exact%3Abtk%20OR%20gene_exact%3Acasp10%20OR%20gene_exact%3Acasp8%29%20AND%20reviewed%3Atrue

def get_UniProt_data(genes):
    # we define the beginning of the url but don't add the genes yet
    # url = 'https://www.uniprot.org/uniprot/?query=organism:%22homo%20sapiens%22%20and%20('
    # url above doesn't work anymore (28.6.2022); use new format instead:
    # url = 'https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=tsv&query=organism_name%3A%22homo%20sapiens%22'
    # url = 'https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=xlsx&query=organism_name%3A%22homo%20sapiens%22'
    url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=tsv&query=organism_name%3A%22homo%20sapiens%22'
    # now we loop over all the genes to add them in our query url
    for gene in genes:
        # for the first gene we use the AND operator, for the rest the OR operator
        # so we add an if else statement like so, to catch the first gene
        if genes.index(gene)==0:
            # gene_url_string = f'gene_exact:{gene.lower()}'
            # url above doesn't work anymore (28.6.2022); use new format instead:
            gene_url_string = f'%20AND%20%28gene_exact%3A{gene.lower()}'
        else:
            # gene_url_string = f'%20or%20gene_exact:{gene.lower()}'
            # url above doesn't work anymore (28.6.2022); use new format instead:
            gene_url_string = f'%20OR%20gene_exact%3A{gene.lower()}' 
        # now we add the part of the url for this gene to the overall url
        url = url + gene_url_string
    # finally we add the other components of our search to the url:
    # search_params = ')%20and%20reviewed:yes&format=tab&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length&sort=score'
    # url above doesn't work anymore (28.6.2022); use new format instead:
    search_params = '%29%20AND%20reviewed%3Atrue&size=500'
    # we combine everything to get our final url (to download data from the uniprot API)
    final_url = url + search_params
    # to download results from uniprot my url should look like this
    # 'https://www.uniprot.org/uniprot/?query=organism:%22homo%20sapiens%22%20and%20(gene_exact:nek1%20or%20gene_exact:kl%20or%20gene_exact:fus)%20and%20reviewed:yes&format=tab&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length&sort=score'
    # apparently this doesn't work anymore (28.6.2022)
    # try this format instead:
    # 'https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=tsv&query=organism_name%3A%22homo%20sapiens%22%20AND%20%28gene_exact%3Abraf%20OR%20gene_exact%3Abrca1%20OR%20gene_exact%3Abrca2%20OR%20gene_exact%3Abtk%20OR%20gene_exact%3Acasp10%20OR%20gene_exact%3Acasp8%29%20AND%20reviewed%3Atrue&size=500'
    # Now we get data from the API and save it as a txt file
    response = requests.get(final_url)
    if response.status_code == 200:
        with open ('uniprot_results.txt', 'w') as file:
            file.write(response.text)
        # now we read this data in as a df
        uniprot_data = pd.read_csv('uniprot_results.txt', sep='\t')
        # and we delete the file uniprot_results.txt again
        os.remove('uniprot_results.txt')
        return uniprot_data
    else:
        # If there is no data, print status code and response
        print(response.status_code, response.text)
        print(f'No data retrieved from UniProt ID for query url {final_url}\n')
        return None # for failure


# Define a function which takes a uniprot ID as input and will download the corresponding mmCif structure from the AlphaFold Database
def get_AlphaFold_structure(uniprot_id):
    mmCIF_url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v2.cif'
    response = requests.get(mmCIF_url)
    if response.status_code == 200:
        data = response.text
        # save the data as mmCIF file
        with open (f'{uniprot_id}_AFprediction.cif', 'w') as f:
            f.write(data)
        print(f'Download of AlphaFold structure for {uniprot_id} successful!')
        return 1 # for success
    else:
        # If there is no data, print status code and response
        print(response.status_code, response.text)
        print(f'No data retrieved for UniProt ID {uniprot_id}\n')
        return 0 # for failure

# ==========================================================================================
# we change to the results directory
os.chdir(results_dir)

# make a folder to store all AlphaFold structures if it doesn't already exist:
alphafold_dir = f'{results_dir}/AlphaFold_structures' #define path
# create Results folder if it doesn't already exist
if not os.path.exists(alphafold_dir):
        os.makedirs(alphafold_dir)

# Now we can use the function above to get data from UniProt for all our gene names and store this data in a dataframe
print('>>> getting UniProt identifiers...')
uniprot_data = get_UniProt_data(genes)
print(f'Complete!\n    Identified {len(uniprot_data)} UniProt enries for a total of {len(genes)} input genes\n')
# We define a df to output at the end of the script, containing relevant information from uniprot and a link to the AlphaFold Database entry
output_df = pd.concat([uniprot_data, pd.DataFrame(columns=['Link AlphaFold Database', 'Status'])])

# We loop over this table and get the Alphafold structure for each UniProt ID (stored in the column 'Entry') using the function defined above
# we will simultaneously fill in the table with additional data (the link to the entry in the AlphaFold database, and the download status)
print('>>> getting AlphaFold2 predicted structures...')
# first we change to the alphafold folder so the structures get saved in the correct place
os.chdir(alphafold_dir)
for index, row in output_df.iterrows():
    uniprot_id = row.Entry
    # download AlphaFold structure
    # the function returns 1 if download was succesful and 0 if it failed (the structure is downloaded in the background)
    # so we use this information to fill in the download status column in the output table
    if get_AlphaFold_structure(uniprot_id) == 1:
        output_df.loc[index, 'Status'] = 'downloaded'
        # we also rename the downloaded structure file (currently in format: uniprotID_AFprediction.cif)
        # it should also contain the gene name (one or more gene names for each structure stored in column output_df['Gene names']
        gene_name = [gene_name for gene_name in output_df.loc[index, 'Gene Names'].split(' ') if gene_name in genes][0]
        os.rename(f'{uniprot_id}_AFprediction.cif', f'{gene_name}_{uniprot_id}_AFprediction.cif')
    else:
        output_df.loc[index, 'Status'] = 'failed'
    # we also add the link to this Entry in the AlphaFold Database to the table
    output_df.loc[index, 'Link AlphaFold Database'] = f'https://www.alphafold.ebi.ac.uk/entry/{uniprot_id}'
print('Complete!\n    All AlphaFold structures have been dowloaded.')
# ==========================================================================================

# write output to a csv file in the results directory
os.chdir(results_dir)
output_df.to_csv('08_AlphaFold_structures.csv', index = False)

# change back to target directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'    o      A total of {len(uniprot_data)} UniProt IDs have been found for the inputted {len(genes)} genes.')
print(f'    o      AlphaFold predicted structures have downloaded for {len(output_df[output_df.Status == "downloaded"])} out of {len(uniprot_data)} identified UniProt IDs.\n')

print('The following files have been created:')
print('   o      08_AlphaFold_structures.csv      (contains information on downloaded AlphaFold structures)\n\n')

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
