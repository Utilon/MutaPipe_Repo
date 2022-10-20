# Parsing clinvar data per gene 
# ****************************
# This script takes no csv file as input, but will automatically read in and parse the xml batch files downloaded from ClinVar
# in the previous script (06_a_download_ClinVar_Data.py). It will:
#      - parse xml variant files (batches of 250) and create a df with ClinVar information for all variants for all input genes
#      - outputs the following files:
#                - a csv file called 06_b_ClinVar_Annotations.csv containing ClinVar Annotations for all variants in all genes
# ===========================================================================================

# Set up
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
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
ap = argparse.ArgumentParser(description="""**** This script takes no csv file as input, but will automatically read in and parse
the xml batch files downloaded from ClinVar in the previous script (06_a_download_ClinVar_Data.py).
It will:
1. loop over all the xml files containing variant information from ClinVar which have been downloaded in script 06_a_download_ClinVar_data.py
2. parse the previously downloaded xml files and create a df with ClinVar information for all variants for all input genes
3. output a csv file called 06_b_ClinVar_Annotations.csv which lists all ClinVar annotations for all variants in all genes       ***""")

ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'specify target directory, default = {target_directory}')
ap.add_argument("-w", "--web_run", type=str2bool, required = False, help=f'Indicate whether MutaPipe is run via a webserver (True) or not (False), default = {str(mutafy_directory)}')
ap.add_argument("-m", "--mutafy", required = False, help=f'set path to mutafy directory where information from previous runs is stored, default = {mutafy_directory}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
web_run = web_run if args["web_run"] == None else args["web_run"]
mutafy_directory = mutafy_directory if args["mutafy"] == None else args["mutafy"]

# ----------------------------------------------------------------------------------------------------------------------------------
# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# ----------------------------------------------------------------------------------------------------------------------------------
#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_06_b.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_06_b.txt', 'a')

# print nice title
print('===============================================================================')
print('*****    Parsing ClinVar Annotations for Input Genes    *****')
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


# change to ClinVar_Annotations  folder (where we have the data of all variants for all genes in xml files; from running the previous script 06_a)
if web_run:
    clinvar_dir = f'{mutafy_directory}/ClinVar_Annotations'
else:
    clinvar_dir = f'{results_dir}/ClinVar_Annotations'
os.chdir(clinvar_dir)

# make a df to populate with relevant info from clinvar xml output for all genes
clinvar_data = pd.DataFrame(columns=['input_gene', 'gene', 'accession', 'title', 'variant_type',  'protein_change', 'aliases', 'clinical significance', 'last_evaluated', 'review_status', 'associated_traits', 'dbs_and_accessions'])

# count how many batches there are per gene and how many genes we have data for:
all_ClinVar_files = [f for f in listdir(clinvar_dir) if isfile(join(clinvar_dir, f))]
# now we have a list of all files. these include a txt file for each gene with the ids for all variants in that gene and xml files batched for 250 variants per file

id_files = [f for f in all_ClinVar_files if 'ids.xml' in f]

# we have data for this number of genes:
n_genes_clinvar_data = len(id_files)

# we create an empty df to populate with all gene names and the corresponding number of data batches (each up to 250 variants)
genes_batches = pd.DataFrame(columns=['gene', 'n_data_batches'])

# we can loop over the list with the id files to parse all corresponding xml batches
for id_file in id_files:
    # the id file is a string in the fomat '06_a_ClinVar_SOD1_ids.txt'
    # we get the gene name by replacing the beginning and the end of the string
    this_gene = id_file.replace('06_a_ClinVar_', '')
    this_gene = this_gene.replace('_ids.xml', '')
    
    # now we get a list of all the corresponding batch xml files for this gene:
    batch_files = [f for f in all_ClinVar_files if f'{this_gene}_data_batch' in f]
    
    # in order to populate the genes_batches df, we get the length of the batch_files list we just created
    genes_batches.loc[len(genes_batches)] = [this_gene, len(batch_files)]
    
    print(f'\n***Gene: {this_gene}***')
    print(f'There are {len(batch_files)} xml batch files for {this_gene}')
    
    # now we loop over the batch files and parse them 
    for batch_file in batch_files:
        print(f'    >>> parsing xml batch {batch_files.index(batch_file)+1} of {len(batch_files)} for {this_gene}')
        # we read it in to extract information
        tree = ET.parse(batch_file)
        root = tree.getroot()
        
        # figured out way to get specific variants and their attributes by trial and error, ahaha
        # we can acces elements of the tree with indexing the root, e.g. root[0][1]
        # [0] first number is alway 0 (DocumentSummarySet ; all records are stored in this first and only set)
        # [0-n] second number a DocumentSummary for a specific record (for a clinvar id)
        # [] third number accesses attributes within a record for one clinvar id
        # so for the the following attributes of the first DocumentSummary, you could do:
        # variant_type = root[0][1][0].text
        # alleles_type = root[0][1][3].text
        # protein_change = root[0][1][15].text

        # but we don't just want to do this for one DocumentSummary, but for all of them
        # so we loop over all DocumentSummaries like this:
        for summary in root.iter('DocumentSummary'):
            sum_gene = summary[9].text
            accession = summary[1].text
            title = summary[3].text
            variant_type = summary[0].text
            protein_change = summary[15].text
            clinical_significance = summary[7][0].text
            last_evaluated = summary[7][1].text
            review_status = summary[7][2].text
            
            all_aliases = []
            for variation_set in summary.iter('variation_set'):
                for variation in variation_set.iter('variation'):
                    for alias in variation.iter('aliases'):
                        aliases = [x.text for x in alias.iter('string')]
                        for a in aliases:
                            all_aliases.append(a)
            # bring aliases in nice format:
            if len(all_aliases)==1:
                all_aliases = all_aliases[0]
            elif len(all_aliases) == 0:
                all_aliases = ''              
                            
            all_associated_traits = []
            dbs_and_accessions = {}
            for trait_set in summary.iter('trait_set'):
                for trait in trait_set.iter('trait'):
                    trait_name = trait[1].text
                    all_associated_traits.append(trait_name)
                    for references in trait.iter('trait_xrefs'):
                        for db in references.iter('trait_xref'):
                            db_name = db[0].text
                            db_accession = db[1].text
                            dbs_and_accessions[db_name] = db_accession
            
            # now we append all these variables to the df clinvar_data
            # however, sometimes the retrieved variants are associated with other genes (not the current query gene)
            # we only append a row to the df if it's for the correct gene (the current query gene)
            if this_gene == sum_gene:
                clinvar_data.loc[len(clinvar_data)] = [this_gene, sum_gene, accession, title, variant_type, protein_change, all_aliases, clinical_significance, last_evaluated, review_status, all_associated_traits, dbs_and_accessions]
                
                # whenever, we create a new row in the clinvar_data df, we write the current version to a csv file.
                # this file will be overwritten everytime a new row is added.
                clinvar_data.to_csv(f'{results_dir}/06_b_ClinVar_Annotations.csv', index = False)
                    
# write df to csv
clinvar_data.to_csv(f'{results_dir}/06_b_ClinVar_Annotations.csv', index = False)

# we also want to write gene specific files to all the gene folders
# we only do this if this is not a webrun
if not web_run:
    # we can loop over all the gene names in the the clinvar_data df
    for gene in clinvar_data.gene.unique():
        # get a slice of the clinvar_data df for only this gene
        clinvar_slice = clinvar_data[clinvar_data.input_gene == gene]
        # we get the name of the gene folder of the current gene like so
        # gene_folder = [name for name in os.listdir(results_dir) if name.startswith(f'{gene}_')][0]
        gene_folder = [name for name in os.listdir(results_dir) if (name.startswith(f'{gene}_')) & ('structures' in name) & ('.csv' not in name)][0]
        # and we write the clinvar slice df to a csv file in the gene folder
        clinvar_slice.to_csv(f'{results_dir}/{gene_folder}/{gene}_06_b_ClinVar_Annotations.csv', index=False)

# Now, if this is a webrun, we will delete all the data downloaded from ClinVar again and only keep the parsed outputs
# contained in the csv files we wrote.
# we do this because this script to parse the ClinVar data will parse all ClinVar data in the set folder
# as there are  continuosly more variants available in ClinVar, we have to download the data everytime anyway,
# so it doesn't make sense to create a directory and add to it for every mutafy run (I think)

# we delete all the data from ClinVar if this is a webrun
if web_run:
    # we change to the Clinvar folder
    os.chdir(clinvar_dir)
    # and we delete all the ClinVar files, they are stored in the list all_ClinVar_files
    # we can do this with a list comprehension
    [os.remove(file) for file in all_ClinVar_files]

# change back to target directory
os.chdir(target_directory)

print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Parsed all ClinVar data for a total of {n_genes_clinvar_data} genes.')

if not web_run:
    print('\nThe following files have been created for each gene and stored in the respective gene folder:')
    print('   o      GENENAME_06_b_ClinVar_Annotations.csv')
    print('            (lists all ClinVar annotations for all variants in this gene)')

print('\nThe following files have been created and stored in the Results folder:')
print('   o      06_b_ClinVar_Annotations.csv')
print('            (lists all ClinVar annotations for all variants in all genes)')

# print script name to console/log file
print(f'\nend of script {script_name}')

# store current date and time in an object and print to console / write to log file
end_time = datetime.now()
print(f'start: {start_time}')
print(f'end: {end_time}\n\n')
print('........................................................................................................................................................\n\n\n')

# close search log
if create_search_log == True:
    sys.stdout.close()