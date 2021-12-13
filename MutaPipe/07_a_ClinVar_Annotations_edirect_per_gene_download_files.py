# Download clinvar data per gene (instead of sending a request for every gene-mismatch)
# ************************************************************************************
# This script takes a csv file (00_search_overview_availability.csv) containing the gene names of all available genes
#  and unavailable genes (ClinVar genes with/without PDB data) as input and will:
#      - download xml files with ids for all variants in ClinVar for each gene
#      - parse xml files and download/retrieve data from ClinVar for all variant ids associated with a gene
#      - parse xml files and create a df with ClinVar information for all variants for all input genes
#      - outputs the following files:
#                - all xml files downloaded from ClinVar are stored in the newly created ClinVar_Annotations folder
#                - a txt file called 07_a_ClinVar_Annotations_genes_no_data_retrieved.txt which lists all genes for which no ClinVar annotations could be retrieved

# Note: Initially, I would also parse the xml files in the same script, however - this took too long on rosalind when
# I wanted to do it with >2500 genes (from ClinVar), therefore I had to split the script.
# Now it's 07a to download information from ClinVar and 07b to parse the downloaded files.
# The code not used in this part of the script has been commented out below
# ===========================================================================================

# Set up
import pandas as pd
import requests
import math
import os
import sys
import argparse
from datetime import datetime

import xml.etree.ElementTree as ET
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
ap = argparse.ArgumentParser(description="""****   This script takes a csv file (00_search_overview_availability.csv) containing the gene names of all available genes and unavailable genes (with/without PDB data) as input and will:
1. download xml files with ids for all variants in ClinVar for each gene
2. parse xml files and download/retrieve data from ClinVar for all variant ids associated with each gene
3. output the following files:
(a) all xml files downloaded from ClinVar are stored in the newly created folder ClinVar_Annotations 
(b) a txt file called 07_a_ClinVar_Annotations_genes_no_data_retrieved.txt which lists all genes for which no ClinVar annotations could be retrieved       ***""")

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
    with open(f'{results_dir}/search_log_00.txt', 'w') as search_log:
        search_log.write(f'Search log for 00_search_pdb.py\n\n')
    sys.stdout = open(f'{results_dir}/search_log_00.txt', 'a')

# store current date and time in an object and print to console / write to log file
start_time = datetime.now()
print(f'start: {start_time}\n')

# ----------------------------------------------------------------------------------------------------------------------------------

# read in data
genes_df = pd.read_csv(f'{results_dir}/00_search_overview_availability.csv')
# filter out genes that have available data
avail_genes = genes_df[genes_df.data_available == True].reset_index()

# make a new folder in current working directory called ClinVar_Annotations
if not os.path.exists(f'{results_dir}/ClinVar_Annotations'):
    os.makedirs(f'{results_dir}/ClinVar_Annotations')
    
# change to ClinVar_Annotations  folder
clinvar_dir = f'{results_dir}/ClinVar_Annotations'
os.chdir(clinvar_dir)

# make a df to populate with relevant info from clinvar xml output for all genes
clinvar_data = pd.DataFrame(columns=['input_gene', 'gene', 'accession', 'title', 'variant_type',  'protein_change', 'aliases', 'clinical significance', 'last_evaluated', 'review_status', 'associated_traits', 'dbs_and_accessions'])


# Set up edirect request to clinvar
# ------------------------------------------
# for any edirect requests
base_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

#create a variable to keep track of genes for which we cannot retrieve data for whatever reason:
genes_no_data_retrieved = []
gene_counter = 0
#  Now we need to initiate the loop over the avail_genes df to get each gene name to query ClinVar:
for gene in avail_genes.gene_name:
    gene_counter += 1
    # with this url I can get all the identifiers of all variants for a gene
    esearch_url = base_url + f'esearch.fcgi?db=clinvar&term={gene}[gene_name]&retmax=100000'
    # esearch_url = base_url + f'esearch.fcgi?db=clinvar&term={gene}[gene_name]' ## we also need to set retmax, otherwise it will only retrieve the first 20 entries or so 
    # url = 'https://www.ncbi.nlm.nih.gov/clinvar/?term=sod1%5Bgene%5D%23' #this doesn't work... unfortunately

    # DOWNLOAD CLINVAR DATA FOR THIS GENE
    # *******************************************
    print(f'\nUsing edirect to query ClinVar for gene {gene_counter} of {len(avail_genes)}: {gene}')

    # try to get ids for this gene from clinvar
    response = requests.get(esearch_url)

    if response.status_code == 200:
        print(f'    Downloading ClinVar ids for variants in {gene}')
        with open(f'07_a_ClinVar_{gene}_ids.xml', 'w') as file:
            file.write(response.text)
    else:
        # If there is no data, print status code and response
        print(response.status_code, response.text)
        print(f'No data file retrieved for {gene}\n')
        # add this gene to the missing genes list
        genes_no_data_retrieved.append(gene)
        continue
   
    # now we can use this file to download the data for these identifiers 
    # so we parse the file to get the ids
    tree = ET.parse(f'07_a_ClinVar_{gene}_ids.xml')
    root = tree.getroot()
        
    # save relevant information from xml file into variables
    identifiers = [x.text for x in root.iter('Id')]
    print(f'    --> There are {len(identifiers)} ClinVar identifiers associated with {gene}')

    # make a string of identiefiers (from the list identifiers) which we can insert into a link to get esummary
    # string must be without commas
    # don't need this code anymore, it's in the paragraph below in the for loop
#     string_identifiers = ''
#     for identifier in identifiers:
#         string_identifiers += (identifier + ',')
#     # remove the last comma
#     string_identifiers = string_identifiers[:-1]

    # problem: sometimes there are so many identifiers, we can get problems with the url getting too long
    # solution: we send only 250 identifiers at once
    # it seems the url cannot be too long, e.g. for AARS1, there are 695 identifiers and if we send all of them in
    # one query string, it doesn't work.
    # I have tried and it works when we send a string with less than 4500 characters
    # identifiers usually have up to 7 characters from what I've seen, so taking 250 at once should be fine (7*250 = 1'750)
    batch_counter = 0
    print(f'    Downloading ClinVar data for {gene} in batches of 250 identifiers per query')
    for i in range(0, len(identifiers), 250):
        batch_counter += 1
        string_identifiers = ''
        for identifier in identifiers[i:i+250]:
            string_identifiers += (identifier + ',')
        # remove the last comma
        string_identifiers = string_identifiers[:-1]            
     
        # now we create a url to get the summary of all our identifiers for this gene
        esummary_url = base_url + f'esummary.fcgi?db=clinvar&id={string_identifiers}'

        # and we get the data from clinvar
        response = requests.get(esummary_url)

        if response.status_code == 200:
            print(f'        Downloading ClinVar data for {gene} variants: Batch {batch_counter} of {math.ceil(len(identifiers)/250)}')
            # new added because sometimes I get an error that the unicode can't be decoded: we write bytes to file instead.
            # previous code is commented out, but not deleted! here:
            # actually, for some reason this takes aaages on rosalind, so instead of just writing bytes as a default, I added
            # a try and except statement to try write the data normally ('w') (previous code) except there is an error,
            # in which case we write bytes ('wb')
            try:
                with open(f'07_a_ClinVar_{gene}_data_batch_{batch_counter}_of_{math.ceil(len(identifiers)/250)}.xml', 'w') as file:
                    file.write(response.text)
            except:
                with open(f'07_a_ClinVar_{gene}_data_batch_{batch_counter}_of_{math.ceil(len(identifiers)/250)}.xml', 'wb') as file:
                    file.write(response.text.encode('utf-8', 'ignore'))
        else:
            # If there is no data, print status code and response
            print(response.status_code, response.text)
            print(f'No data file retrieved for {gene}: Batch {batch_counter} of {math.ceil(len(identifiers)/250)}\n')
            # add this gene to the missing genes list
            genes_no_data_retrieved.append(gene)
            continue
            
#         # yippiee, now we have the data of all variants for this gene in an xml file!
#         # we read it in to extract information
#         tree = ET.parse(f'07_a_ClinVar_{gene}_data_batch_{batch_counter}_of_{math.ceil(len(identifiers)/250)}.xml')
#         root = tree.getroot()
# 
#         # figured out way to get specific variants and their attributes by trial and error, ahaha
#         # we can acces elements of the tree with indexing the root, e.g. root[0][1]
#         # [0] first number is alway 0 (DocumentSummarySet ; all records are stored in this first and only set)
#         # [0-n] second number a DocumentSummary for a specific record (for a clinvar id)
#         # [] third number accesses attributes within a record for one clinvar id
#         # so for the the following attributes of the first DocumentSummary, you could do:
#         # variant_type = root[0][1][0].text
#         # alleles_type = root[0][1][3].text
#         # protein_change = root[0][1][15].text
# 
#         # but we don't just want to do this for one DocumentSummary, but for all of them
#         # so we loop over all DocumentSummaries like this:
#         for summary in root.iter('DocumentSummary'):
#             sum_gene = summary[9].text
#             accession = summary[1].text
#             title = summary[3].text
#             variant_type = summary[0].text
#             protein_change = summary[15].text
#             clinical_significance = summary[7][0].text
#             last_evaluated = summary[7][1].text
#             review_status = summary[7][2].text
#             
#             all_aliases = []
#             for variation_set in summary.iter('variation_set'):
#                 for variation in variation_set.iter('variation'):
#                     for alias in variation.iter('aliases'):
#                         aliases = [x.text for x in alias.iter('string')]
#                         for a in aliases:
#                             all_aliases.append(a)
#             # bring aliases in nice format:
#             if len(all_aliases)==1:
#                 all_aliases = all_aliases[0]
#             elif len(all_aliases) == 0:
#                 all_aliases = ''              
#                             
#             all_associated_traits = []
#             dbs_and_accessions = {}
#             for trait_set in summary.iter('trait_set'):
#                 for trait in trait_set.iter('trait'):
#                     trait_name = trait[1].text
#                     all_associated_traits.append(trait_name)
#                     for references in trait.iter('trait_xrefs'):
#                         for db in references.iter('trait_xref'):
#                             db_name = db[0].text
#                             db_accession = db[1].text
#                             dbs_and_accessions[db_name] = db_accession
#             
#             # now we append all these variables to the df clinvar_data
#             # however, sometimes the retrieved variants are associated with other genes (not the current query gene)
#             # we only append a row to the df if it's for the correct gene (the current query gene)
#             if gene == sum_gene:
#                 clinvar_data.loc[len(clinvar_data)] = [gene, sum_gene, accession, title, variant_type, protein_change, all_aliases, clinical_significance, last_evaluated, review_status, all_associated_traits, dbs_and_accessions]
#                 
#                 # whenever, we create a new row in the clinvar_data df, we write the current version to a csv file.
#                 # this file will be overwritten everytime a new row is added.
#                 clinvar_data.to_csv(f'{results_dir}/07_a_ClinVar_Annotations.csv', index = False)
#                 # we also write a textfile containing all gene names which were not available for ClinVar
#                 with open(f'{results_dir}/07_a_ClinVar_Annotations_genes_no_data_retrieved.txt', 'w') as file:
#                     for element in genes_no_data_retrieved:
#                         file.write(element + '\n')
#                     
# 
# # write df to csv
# clinvar_data.to_csv(f'{results_dir}/07_a_ClinVar_Annotations.csv', index = False)
# 
# # now we can use the columns protein_change and aliases to see if the mismatches from the pdb are in clinvar!!!!!
# # then we can directly get the information we need for each mismatch from the clinvar_data df without having to
# # query clinvar for each individual mismatch! :-) --> in separate script!

# we also write a textfile containing all gene names which were not available for ClinVar
with open(f'{results_dir}/07_a_ClinVar_Annotations_genes_no_data_retrieved.txt', 'w') as file:
    for element in genes_no_data_retrieved:
        file.write(element + '\n')


# change back to target directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Queried ClinVar for a total of {len(avail_genes)} genes using edirect.')

print('\n     All xml outputs from ClinVar stored in .xml format in the folder ClinVar_Annotations\n')

print('The following files have been created and stored in the Results folder:')
print('   o      07_a_ClinVar_Annotations_genes_no_data_retrieved.txt       (lists all genes for which no ClinVar annotations could be retrieved)\n\n')


# store current date and time in an object and print to console / write to log file
end_time = datetime.now()
print(f'start: {start_time}\n')
print(f'end: {end_time}\n\n')

# close search log
if create_search_log == True:
    sys.stdout.close()
