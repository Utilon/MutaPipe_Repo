# Download clinvar data for all genes with available PDB data
# ************************************************************************************
# This script takes a csv file (00_search_overview_availability.csv) containing the gene names of all available genes
#  and unavailable genes (ClinVar genes with/without PDB data) as input and will:
#      - download xml files with ids for all variants in ClinVar for each available gene
#      - parse xml files and download/retrieve data from ClinVar for all variant ids associated with a gene
#      - outputs the following files:
#                - all xml files containing variant information downloaded from ClinVar are stored in the newly created ClinVar_Annotations folder
#                - a txt file called 06_a_ClinVar_Annotations_genes_no_data_retrieved.txt which lists all genes for which no ClinVar annotations could be retrieved
# ===========================================================================================

# Set up
import pandas as pd
import requests
import math
import os
from os.path import exists
import sys
import argparse
from datetime import datetime

import xml.etree.ElementTree as ET

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
web_run = True # specify if pdb mmcif and fasta files should be stored in separate directory
mutafy_directory = f'{target_directory}/mutafy' # set path to folder where structures will be/are stored        
                                                                                      
# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****   This script takes a csv file (00_search_overview_availability.csv) containing the gene names of all available genes and unavailable genes (with/without PDB data) as input and will:
1. download xml files with ids for all variants in ClinVar for each gene
2. parse xml files and download/retrieve data from ClinVar for all variant ids associated with each gene
3. output the following files:
(a) all xml files downloaded from ClinVar are stored in the newly created folder ClinVar_Annotations 
(b) a txt file called 06_a_ClinVar_Annotations_genes_no_data_retrieved.txt which lists all genes for which no ClinVar annotations could be retrieved       ***""")

ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')
ap.add_argument("-w", "--web_run", type=str2bool, required = False, help=f'Indicate whether MutaPipe is run via a webserver (True) or not (False), default = {str(mutafy_directory)}')
ap.add_argument("-m", "--mutafy", required = False, help=f'set path to mutafy directory where information from previous runs is stored, default = {mutafy_directory}')

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
    with open(f'{results_dir}/search_log_06_a.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_06_a.txt', 'a')

# print nice title
print('===============================================================================')
print('*****    Downloading ClinVar Annotations for Input Genes    *****')
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

# read in data
genes_df = pd.read_csv(f'{results_dir}/00_search_overview_availability.csv')
# filter out genes that have available PDB data
avail_genes = genes_df[genes_df.data_available == True].reset_index()

# CHECK FOR WHICH GENES WE NEED TO DOWNLOAD CLINVAR DATA
# read in mutafy data (06_b_ClinVar_Annotations_mutafy.csv) if available
# we will use this to see if ClinVar data has already been downloaded for the genes of interest in a previous webrun
if web_run:
    if exists(f'{mutafy_directory}/06_b_ClinVar_Annotations_mutafy.csv'):
        annotations_mutafy = pd.read_csv(f'{mutafy_directory}/06_b_ClinVar_Annotations_mutafy.csv')
        # for now, we only download the ClinVar annotations once per gene (newly published variants will be neglegted)
        # so we will use the annotations_mutafy to see if we have data for any of the avail_genes already
        # and we update avail_genes accordingly (so it only contains new genes for which ClinVar data still hast to be downloaded)
        # get list of unique genes for which data is already available
        genes_with_ClinVar_data = list(annotations_mutafy.input_gene.unique())
        # get list of unique genes of current run
        genes_current_run = list(avail_genes.gene_name.unique())
        # now we update the list avail_genes, which stores all the genes for which ClinVar data will be downloaded in the
        # next step. For this, we drop all the genes in the genes_current_run list if they are already in the genes_with_ClinVar_data
        # list
        new_genes = [gene for gene in genes_current_run if gene not in genes_with_ClinVar_data]
        avail_genes = avail_genes[avail_genes.gene_name.isin(new_genes)]


# make a new folder directory called ClinVar_Annotations
if web_run:
    clinvar_dir = f'{mutafy_directory}/ClinVar_Annotations'
else:
    clinvar_dir = f'{results_dir}/ClinVar_Annotations'
    
if not os.path.exists(clinvar_dir):
    os.makedirs(clinvar_dir)
    
# change to ClinVar_Annotations  folder
os.chdir(clinvar_dir)

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
        with open(f'06_a_ClinVar_{gene}_ids.xml', 'w') as file:
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
    tree = ET.parse(f'06_a_ClinVar_{gene}_ids.xml')
    root = tree.getroot()
        
    # save relevant information from xml file into variables
    identifiers = [x.text for x in root.iter('Id')]
    print(f'    --> There are {len(identifiers)} ClinVar identifiers associated with {gene}')

    # make a string of identiefiers (from the list identifiers) which we can insert into a link to get esummary
    # string must be without commas
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
                with open(f'06_a_ClinVar_{gene}_data_batch_{batch_counter}_of_{math.ceil(len(identifiers)/250)}.xml', 'w') as file:
                    file.write(response.text)
            except:
                with open(f'06_a_ClinVar_{gene}_data_batch_{batch_counter}_of_{math.ceil(len(identifiers)/250)}.xml', 'wb') as file:
                    file.write(response.text.encode('utf-8', 'ignore'))
        else:
            # If there is no data, print status code and response
            print(response.status_code, response.text)
            print(f'No data file retrieved for {gene}: Batch {batch_counter} of {math.ceil(len(identifiers)/250)}\n')
            # add this gene to the missing genes list
            genes_no_data_retrieved.append(gene)
            continue
            
# we also write a textfile containing all gene names which were not available for ClinVar
with open(f'{results_dir}/06_a_ClinVar_Annotations_genes_no_data_retrieved.txt', 'w') as file:
    for element in genes_no_data_retrieved:
        file.write(element + '\n')

# change back to target directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Queried ClinVar for a total of {len(avail_genes)} genes using edirect.')

print('\n     All xml outputs from ClinVar stored in .xml format in the folder ClinVar_Annotations\n')

print('The following files have been created and stored in the Results folder:')
print('   o      06_a_ClinVar_Annotations_genes_no_data_retrieved.txt')
print('            (lists all genes for which no ClinVar annotations could be retrieved)\n\n')

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