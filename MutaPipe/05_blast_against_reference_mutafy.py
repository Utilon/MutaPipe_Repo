# This script takes a csv file (04_fasta_combined_info.csv) containing information for all genes extracted from their
#  respective fasta files and fasta_ex files as input and will:
#      - extract the reference sequence for the gene from the uniprot fasta file
#      - perform BLASTp with each sequence of each structure (=each row of the df) against the respective
#        reference sequence (e.g. FUS canonical sequence for all sequences in all FUS structures)
#      - uses the blast output (xml files) to identify mismatches and add this information to the df
#      - outputs the following files:
#                - all reference sequences and unique sequences of all structures are stored in .fasta format in the newly created Results/RefSeqs folder
#                - all BLASTp outputs are stored in .xml format in the newly created Results/RefSeqs folder
#                - a csv file called 05_blastp_results.csv listing all the information in the input file and the corresponding blastp results
# ===========================================================================================

# Set up
import pandas as pd
import numpy as np
import os
from os import listdir
from os.path import isfile, join, exists
import ast
import re
import sys
import argparse
from datetime import datetime

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
blastp_path =  '/usr/local/ncbi/blast/bin/blastp' # First run 'which blastp' on CommandLine to find the full path
                                                                         #  to blastp and give that as argument to NcbiblastpCommandline.
target_directory = os.getcwd()    # set target directory (where Results folder is located)

# figure out path to uniprot fasta (no argument) to set the default value for uniprot_fasta (see one line further down)
MutaPipe_dir = os.path.dirname(os.path.realpath(__file__))
uniprot_fasta = f'{MutaPipe_dir}/Uniprot_reference_seqs/UP000005640_9606.fasta' # specify default path to uniprot reference fasta

web_run = True # specify if pdb mmcif and fasta files should be stored in separate directory
mutafy_directory = f'{target_directory}/mutafy' # set path to folder where structures will be/are stored        

# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****    This script takes a csv file (04_fasta_combined_info.csv) containing information for all input genes extracted from their
#  respective fasta files and fasta_ex files as input and will: 
1. retrieve the reference sequence for the gene associated with each structure from the uniprot reference fasta file 
2. perform BLASTp with each sequence of each structure against the respective reference sequence (e.g. FUS canonical sequence for all sequences in all FUS structures) 
3. uses the blast output (xml files) to identify mismatches / variants in all sequences 
4. outputs the following files: 
(1) all reference sequences and unique sequences of all structures are stored in .fasta format in the newly created Results/RefSeqs folder 
(2) all BLASTp outputs are stored in .xml format in the newly created Results/RefSeqs folder 
(3) a csv file called 05_blastp_results.csv listing all the information on all sequences in all structures for all genes and the corresponding blastp results    ***""")

ap.add_argument('-bp','--blastp_path', required=False, help=f"Specify the path to blastp on your system ; default = {blastp_path}")
ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'Write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'Specify target directory, default = {target_directory}')
ap.add_argument("-refseq", "--reference_sequences", required = False, help=f'Specify path to uniprot reference fasta, default = {uniprot_fasta}')
ap.add_argument("-w", "--web_run", type=str2bool, required = False, help=f'Indicate whether MutaPipe is run via a webserver (True) or not (False), default = {str(mutafy_directory)}')
ap.add_argument("-m", "--mutafy", required = False, help=f'set path to mutafy directory where information from previous runs is stored, default = {mutafy_directory}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
blastp_path = blastp_path if args["blastp_path"] == None else args["blastp_path"]
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
uniprot_fasta = uniprot_fasta if args["reference_sequences"] == None else args["reference_sequences"]
web_run = web_run if args["web_run"] == None else args["web_run"]
mutafy_directory = mutafy_directory if args["mutafy"] == None else args["mutafy"]
# ----------------------------------------------------------------------------------------------------------------------------------
# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# ----------------------------------------------------------------------------------------------------------------------------------
#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_05.txt', 'w') as search_log:
        search_log.write(f'Search log for {script_name}\n\n')
    sys.stdout = open(f'{results_dir}/search_log_05.txt', 'a')
    
# print nice title
print('===============================================================================')
print('*****    BLASTp against Reference Sequence for all Identified PDB Sequences   *****')
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

# we make a new folder in the Results folder to store all reference sequences (permanently) and individual fasta files (temporarily)
# which we need during our for loop when doing the blastp
#  in case this is a webrun, the RefSeqs folder will be in the mutafy_directory
if web_run:
    refseqs_dir = f'{mutafy_directory}/RefSeqs'
else:
    refseqs_dir = f'{results_dir}/RefSeqs'
if not os.path.isdir(refseqs_dir):
    os.mkdir(f'{refseqs_dir}')
    
# change to refseqs_dir
os.chdir(refseqs_dir)

# Read in data from csv files
fasta_df = pd.read_csv(f'{results_dir}/04_fasta_combined_info.csv')

# if this is a web run, we also read in the csv file from script 01 which lists all the new pdb ID's to be downloaded / parsed / blasted!
if web_run:
    df_new_structures_to_blast = pd.read_csv(f'{mutafy_directory}/01_new structures_to_be parsed_mutafy.csv')

# create an empty list to populate with results from blastp
xml_list = []
# we also create empty lists to capture all our WARNINGS (as defined by me, haha) and write them to a csv file at the end
refseq_warnings_list = []
blastp_warnings_list = []

# we also read in the fasta file downloaded from uniprot containing all canonical sequences for the human genome:
uniprot_refseqs = list(SeqIO.parse(uniprot_fasta, 'fasta'))
# the gene_name is included in the description of the fasta record with the canonical sequence: 
# e.g. 'GN=SOD1'  for SOD1

# define variable to keep track of previous gene and previous reference sequence and the name of the reference sequence fasta file,
# so we only have to identify the reference sequence once per gene (and not for every row in the fasta_df)
previous_gene = 'no_gene'
previous_refseq = None
previous_refseq_fasta = 'no_fasta'
gene_counter = 0

# now we want to loop over the fasta_df and for each row perform a blast p against the reference sequence for that gene
for fasta_df_index, row in fasta_df.iterrows():
    gene = row.gene_name
    structure_id = row.structure_id
    
    # if this is a webrun, we first check if there are any new fasta files which have to be blasted (haha) in this folder/for this gene:
    if web_run:
            # we check if the gene and structure id for the current loop is found in the mutafy df of new structures
            # if this is not the case, it means this is not a new structure, so we continue / skip the rest of the loop for this row
            # as the BLASTp has already been performed in a previous mutafy run for this structure
            # first we check if the gene is in the mutafy df listing new structures
            if not (gene in df_new_structures_to_blast.gene.to_list()):
                continue
            # and then we check if the the specific structure id is in the mutafy df listing new structures 
            elif not (structure_id in ast.literal_eval(df_new_structures_to_blast[df_new_structures_to_blast.gene == gene].new_pdb_ids.values[0])):
                continue            

    sequence = Seq(row.sequence)
    # in order to store this sequence properly, we need to create a sequence record, see: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec35
    # Including an identifier is very important if we want to output the SeqRecord to a file.
    seq_rec = SeqRecord(sequence, id=f'{gene}_{row.structure_id}_{row.chain_name.replace(" ", "")}')
    
    # in order to perform a blastp with this sequence against the reference sequence,
    # we need both sequences as fasta files, so we make a fasta file of the current sequence record
    # we create a variable containing the name of the fasta file we are going to write later
    # we need to limit this filename to a max number of characters, otherwise, we get an error when trying to write the fasta file
    # I know that up to 76 characters certainly work, so we only take the first 76 characters of the filename in case it is longer than 76 characters
    fasta = f"{gene}_{row.structure_id}_{row.chain_name.replace(' ', '')}"[:76] + ".fasta"
    # we write the fasta file assigning it the filename we just stored in the variable fasta
    SeqIO.write(seq_rec, fasta, 'fasta')
    
    # we get the reference sequence for this gene if we don't already have it:
    if gene != previous_gene:
        gene_counter += 1
        # we get the reference sequence from the uniprot fasta file of canonical sequences by searching for the gene name
        # we make a list and include only sequences/fasta records that contain the gene name in the description
        # in order to get the correct sequence, e.g. NEK1 instead of NEK10, we include a space after the gene name!!!
        potential_reference_sequences = [sequence for sequence in uniprot_refseqs if f'GN={gene + " "}' in sequence.description]
        # there should only be one reference sequence, so we take the first and only element of this list as our reference sequence
        # print WARNING if more than one reference sequence has been identified (in this case the first one is taken as a reference sequence for the next step)
        if len(potential_reference_sequences) > 1:
            print(f'WARNING: more than 1 reference sequence identified for gene {gene}: {len(potential_reference_sequences)} potential reference sequences \n               First identified sequence with GN={gene} will be used to continue for gene {gene}')
            refseq_warnings_list.append({'gene': gene, 'n_refseqs_identified': len(potential_reference_sequences)})
        # but e.g. for SMN1 there's no reference at all (its name is snm2 for some reason!? download new uniprot fasta!)
        # so we add a try and and except statement:
        try:
            reference_sequence = potential_reference_sequences[0]
        except IndexError:
            print(f'WARNING!        No reference sequence found for gene {gene}')
            # if we don't find a reference sequence for this gene, we add info to refseq_warnings_list
            refseq_warnings_list.append({'gene': gene, 'n_refseqs_identified': 0})
            # and reset the variables for previous_refseq, previous_gene, and previous_refseq_fasta
            previous_refseq = None
            previous_gene = gene
            previous_refseq_fasta = 'no_fasta'
            # and continue (skip the rest of the loop for this row)
            continue
        # we write a fasta file for the current reference sequence that we can later use to blastp
        # we define a variable for the name of this file, so we can read it in again later
        refseq_fasta = f'{gene}_reference.fasta'
        # now we write the fasta file
        print(f'>>> Writing canonical reference sequence for {gene} (gene {gene_counter} of {len(fasta_df.gene_name.unique())}) to fasta file for blastp')
        SeqIO.write(reference_sequence, refseq_fasta, 'fasta')
        # we update the variables previous_refseq_fasta, previous_refseq, and previous_gene:
        previous_gene = gene
        previous_refseq = reference_sequence
        previous_refseq_fasta = refseq_fasta
        
    # if this gene is the same as the previous one, the current reference sequence is the same as the previous one too
    elif gene == previous_gene:
        # the reference sequence is the same as the previous one, so we don't need to write a new fasta file (it already exists)
        # in order to prevent potential bugs, we still update our variable reference_sequence for this gene,
        # even if this is potentially redundant as we don't need this variable but the fasta file
        reference_sequence = previous_refseq
        # the reference fasta file is the same as the previous one
        refseq_fasta = previous_refseq_fasta
        # if there is no reference sequence for this gene (because it's not available / couldn't be extracted),
        # we use continue to skip the rest of the loop
        if reference_sequence is None:
            continue
        
    # so now that we have the sequence and the reference sequence in fasta format,
    # we want to perform a blastp of the two
    # useful info on how to do this here: https://www.biostars.org/p/42687/
    # also here: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec127
    
    # =======================================================================    
    
    # Run BLAST and save output as xml file
    print(f'        >>> perfoming blastp for {gene} fasta file {fasta}')
    blastp_cline = NcbiblastpCommandline(cmd=blastp_path, query=fasta, subject= refseq_fasta, outfmt=5, out=fasta.replace('fasta', 'xml'))
    try:
        stdout, stderr = blastp_cline()
    except:
        print(f'        blastp UNSUCCESFUL for {gene} fasta {fasta}')
        # If BLASTP failed, we add a warning to the warning list and continue (skip rest of the loop)
        blastp_warnings_list.append({'gene': gene, 'blast_status': f'FAILED for {fasta}', 'n_alignments': np.nan})
        continue 

    # parse xml file/blastp output to read blastp results
    print(f'        >>> parsing blastp output for fasta file {fasta}\n')
    result_handle = open (fasta.replace('fasta', 'xml'))
    blast_record = NCBIXML.read(result_handle)

    # =======================================================================
    # now we can use the blast_record to extract information on mismatches, gaps etc.
    
    # there is ideally only one alignment in the blast_record (only 1 subject sequence given, i.e. our reference sequence for each gene)
    # or there could be no alignments, so we add a try and except statement
    try:
        alignment = blast_record.alignments[0]
    except IndexError:
        blastp_warnings_list.append({'gene': gene, 'blast_status': f'No alignments for {fasta}', 'n_alignments': 0})
        continue
    
    # in case there is more than one alignment: print warning
    if len(blast_record.alignments) > 1:
        print(f'WARNING: more than one alignment found in blastp output {fasta.replace("fasta", "xml")}')
        print('            NOTE: only first top alignment used in further analysis')
        blastp_warnings_list.append({'gene': gene, 'blast_status': f'>1 alignment for {fasta}', 'n_alignments': len(blast_record.alignments)})
    
    # Get the best high scoring segment pair for the current alignment
    # there may be multiple high scoring segment pairs (HSPs) in an alignment,
    # if that's the case, we only take the one with the highest bit score
    # (otherwise, we'd have to create a new row for each hsp in the output df or add more columns for more hsp data or something)
    hsp = None
    hsp_counter = 0
    best_bit_score = 0
    best_hsp_counter = 0 
    for this_hsp in alignment.hsps:
        hsp_counter += 1
        if this_hsp.bits > best_bit_score:
            best_bit_score = this_hsp.bits
            hsp = this_hsp
            best_hsp_counter = hsp_counter
            
    # =======================================================================
                  
    # find the number of mismatches and close mismatches using the hsp.match sequence
    mismatches = len([letter for letter in hsp.match if letter == ' '])
    close_mismatches = len([letter for letter in hsp.match if letter == '+'])

    # find the positions of mismatches and close mismatches in the hsp.match sequence
    mismatch_pos = [char for char, ltr in enumerate(hsp.match) if ltr == ' ']
    close_mismatch_pos = [char for char, ltr in enumerate(hsp.match) if ltr == '+']

    # define two empty lists to populate with all mismatch and close mismatch AA substitutions                                                   
    mismatch_substitutions = []
    close_mismatch_substitutions = []
            
    # Mismatch substitions
    # now we loop over all mistmach positions to extract the corresponding amino acids
    # from the hsp.query and hsp.sbjct sequences
    for pos in mismatch_pos:
        query_mis = hsp.query[pos]
        sbjct_mis = hsp.sbjct[pos]
        # when counting residues, we start at 1 not at 0 like in python, thus we add +1 to pos
        substitution = sbjct_mis + str(pos+1) + query_mis
        # now we can append this new AA substitution to the mismatch_substitutions list
        # we are not interested in indels at this point, so we only append AA substitutions to the list
        if '-' not in substitution:
            mismatch_substitutions.append(substitution)
            
    # we do the same for close mismatches                        
    for pos in close_mismatch_pos:
        query_mis = hsp.query[pos]
        sbjct_mis = hsp.sbjct[pos]
        substitution = sbjct_mis + str(pos+1) + query_mis
        close_mismatch_substitutions.append(substitution)
                
    # get the number of mismatches excl. gaps and the number of gaps
    mismatches_excl_gaps = len(mismatch_substitutions)
    gaps = hsp.gaps

    # =======================================================================
    # GAPS / INDELS
    
    # gaps and indels in query
    # to get the indices and the length of all the gaps in the query sequence, we do:
    # get indices of all '-' in query sequence
    gap_indices_query = [letter for letter, ltr in enumerate(hsp.query) if ltr =='-']
    # find all continuous gaps/indels in query sequence
    gaps_in_seq_query = re.findall('-*', hsp.query)
    gaps_list_query = [entry for entry in gaps_in_seq_query if '-' in entry]
    # number of indels in query sequence is the length of the list we just created
    indels_query = len(gaps_list_query)
    # now we loop over gaps query list and append the length of each gap
    # to the variable all_gaps_query
    all_gaps_query = 0
    for gap in gaps_list_query:
        all_gaps_query += len(gap)
        
    # in order to get all the corresponding amino acids from the other sequence
    # if gaps have been introduced in one sequence, we use the above defined variables
    # and extract the sequences of interest in a for loop as follows
    
    # first we define an empty list to store all gap opening indices
    gaps_start_pos_query = []
    # start/opening positions are more than 1 integer bigger than the previous gap position
    # we thus set the previous_index_query to -2 (not -1 as python starts counting with index 0)
    previous_index_query = -2
    # create empty dict to populate with {query_seq : sbjct_seq, ...}
    gap_seqs_dict_query = {}
    # for now we set query_seq and sbjct_seq to ''.
    # We will append letters/AAs/gaps here in the for loop below
    query_seq = ''
    sbjct_seq = ''
    # now we loop over all the gap indices (list we created above)
    for index in gap_indices_query:
        # start new gap if more than 1 increase from previous index
        if index - previous_index_query > 1:
            # append previously constructed query_seq and sbjct_seq to dict unless empty
            if query_seq != '':
                gap_seqs_dict_query[query_seq] = sbjct_seq
            # append current index to list of gap opening positions
            gaps_start_pos_query.append(index)
            # find the letters/AAs/gaps in both sequences at this index
            # we set query_seq/sbjct_seq equal to the letter at this index in order to start
            # a new sequence (if it's a longer gap, more letters will be appended, see elif statement)
            query_seq = hsp.query[index]
            sbjct_seq = hsp.sbjct[index]
        # if this is a continuous gap, then the current index is only 1 integer bigger than the previous one
        elif index - previous_index_query == 1:
            # in this case, we do not start a new sequence, but append the letters at this index
            # in the hsps to the query_seq and the sbjct_seq
            query_seq += hsp.query[index]
            sbjct_seq += hsp.sbjct[index]
        # we update the variable previous_index_query
        previous_index_query = index
    # append the last entry to the dict (we do this outside the for-loop because in it,
    # the entry only get appended to dict if the next gap is read, so the last on doesn't get automatically attached)               
    gap_seqs_dict_query[query_seq] = sbjct_seq
    # we also get a list of the lengths of all gaps in the query sequence:
    gaps_len_query = [len(gap) for gap in gaps_list_query]                
    
    # now we do the same as above in the subject sequence!
    # gaps and indels in sbjct
    gap_indices_sbjct = [letter for letter, ltr in enumerate(hsp.sbjct) if ltr =='-']
    gaps_in_seq_sbjct = re.findall('-*', hsp.sbjct)
    gaps_list_sbjct = [entry for entry in gaps_in_seq_sbjct if '-' in entry]
    indels_sbjct = len(gaps_list_sbjct)
    all_gaps_sbjct = 0
    for gap in gaps_list_sbjct:
        all_gaps_sbjct += len(gap)
    
    gaps_start_pos_sbjct = []
    previous_index_sbjct = -2
    # create empty dict to populate with {sbjct_seq : query_seq, ...}
    gap_seqs_dict_sbjct = {}
    query_seq = ''
    sbjct_seq = ''
    for index in gap_indices_sbjct:
        # start new gap if more than 1 increase from previous index
        if index - previous_index_sbjct > 1:
            # append previously constructed query_seq and sbjct_seq to dict unless empty
            if query_seq != '':
                gap_seqs_dict_sbjct[sbjct_seq] = query_seq                        
            gaps_start_pos_sbjct.append(index)
            query_seq = hsp.query[index]
            sbjct_seq = hsp.sbjct[index]
        elif index - previous_index_sbjct == 1:
            query_seq += hsp.query[index]
            sbjct_seq += hsp.sbjct[index]
        previous_index_sbjct = index
        
    # append the last entry to the dict (we do this outside the for-loop because in it, the entry only get appended to dict if the next gap is read, so the last on doesn't get automatically attached in the foor-loop)               
    gap_seqs_dict_sbjct[sbjct_seq] = query_seq
    gaps_len_sbjct = [len(gap) for gap in gaps_list_sbjct]
    
    # now we have all the information from this blastp output stored in variables, we add it to our xml_list:
    xml_list.append({'gene_name': row.gene_name , 'structure_id': row.structure_id, 'chain_name': row.chain_name,
                     'uniprot_id': row.uniprot_id, 'description': row.description, 'species': row.species,
                     'description_ex': row.description_ex, 'sequence': row.sequence, 'alignment_length': alignment.length,
                     'hsp_number': best_hsp_counter, 'hsp_length': hsp.align_length, 'score': hsp.score, 'bit': hsp.bits,
                     'e-value': hsp.expect, 'similarity': hsp.identities / hsp.align_length,
                     'mismatches_incl_gaps': mismatches, 'mismatches_excl_gaps': mismatches_excl_gaps,
                     'close_mismatches': close_mismatches, 'mismatch_substitutions': str(mismatch_substitutions),
                     'close_mismatch_substitutions': str(close_mismatch_substitutions),
                     'gaps': hsp.gaps, 'gaps_in_query': all_gaps_query, 'gaps_in_sbct': all_gaps_sbjct,
                     'gaps_start_pos_query': str(gaps_start_pos_query), 'gaps_length_query': str(gaps_len_query),
                     'query_del_sbcjt_ins': str(gap_seqs_dict_query), 'gaps_start_pos_sbjct': str(gaps_start_pos_sbjct),
                     'gaps_length_sbjct': str(gaps_len_sbjct), 'sbjct_del_query_ins': str(gap_seqs_dict_sbjct),
                     'query_sequence': hsp.query, 'sbjct_sequence': hsp.sbjct, 'match_sequence': hsp.match})
                     
# Now that we have all the information from all the records for all genes, we convert the xml_list
# to a df and call it xml_df, and the warning lists to dfs as well
xml_df = pd.DataFrame(xml_list)
refseq_warnings = pd.DataFrame(refseq_warnings_list)
blastp_warnings = pd.DataFrame(blastp_warnings_list)

# Clean up the RefSeqs folder
# we make a subdirectory in this folder called PDB_seqs_and_blastp_outputs
# we will move all the fasta files from the PDB structures and the blastp outputs there and
# only keep the reference fasta directly in the RefSeqs directory
# first we get a list of all files in the current folder (RefSeqs)
all_fasta_files_for_blastp = [f for f in listdir(refseqs_dir) if isfile(join(refseqs_dir, f))]
# now we get a list of all the files that are not reference fasta files and move them to a new subfolder
PDB_fasta_files_for_blastp = [f for f in all_fasta_files_for_blastp if not 'reference.fasta' in f]
# make the folder
if not os.path.isdir('PDB_seqs_and_blastp_outputs'):
    os.mkdir('PDB_seqs_and_blastp_outputs')
# and move the files
for PDB_fasta in PDB_fasta_files_for_blastp:
    os.rename(f'{PDB_fasta}', f'PDB_seqs_and_blastp_outputs/{PDB_fasta}')              

# now, if this is a webrun, we want to read in the data from previous mutafy runs and
# get a slice for all the genes of the current seach (otherwise only newly blasted structures will be listed in the output file)
# also, we want to update already existing mutafy data and add new structures to it!
if web_run:
    if exists(f'{mutafy_directory}/05_blastp_results_mutafy.csv'):
        mutafy_blastp_results = pd.read_csv(f'{mutafy_directory}/05_blastp_results_mutafy.csv')
        # we drop all the rows from our xml_df for which we didn't try to perform a blastp
        # we take only the rows with the new structure ids (for which we had to perform a blastp)
        # the relevant structure ids are stored in the df_new_structures_to_blast
        ids_to_add = df_new_structures_to_blast.new_pdb_ids.to_list()
        ids_to_add = [ast.literal_eval(e) for e in ids_to_add]
        # now we have a list of list, to flatten it we do the following:
        # see this post: https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
        ids_to_add = [item for sublist in ids_to_add for item in sublist]
        # now we get a slice of the xml_df containing only the data for the ids_to_add (the ones with new blastp results)
        # we also convert all the values to strings for concatenation (next step).
        xml_slice_to_add = xml_df[xml_df.structure_id.isin(ids_to_add)].astype('str')
        # no we can add the xml slice (containig the new blastp results) to the mutafy data from previous runs
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        updated_mutafy_blastp_results = pd.concat([mutafy_blastp_results, xml_slice_to_add], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name and then structure id
        updated_mutafy_blastp_results.sort_values(by=['gene_name', 'structure_id', 'chain_name'], inplace=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_blastp_results.to_csv(f'{mutafy_directory}/05_blastp_results_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the data for all the structures for all genes of the current webrun
        xml_df = updated_mutafy_blastp_results[updated_mutafy_blastp_results.gene_name.isin(list(fasta_df.gene_name.unique()))]
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        xml_df.to_csv(f'{mutafy_directory}/05_blastp_results_mutafy.csv', index = False)
        # we need a variable called xml_slice_to_add here anyway, to get the lenght of it later in a print statement
        # and print the number of sequences blasted in this cycle
        xml_slice_to_add = xml_df
        
    if exists(f'{mutafy_directory}/05_refseq_warnings_mutafy.csv'):
        mutafy_refseq_warnings = pd.read_csv(f'{mutafy_directory}/05_refseq_warnings_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        # in order to do that, we convert all the values in the refseq_warnings df to strings
        refseq_warnings = refseq_warnings.astype('str')
        # now we can concatenate the two dfs
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        updated_mutafy_refseq_warnings = pd.concat([mutafy_refseq_warnings, refseq_warnings], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name 
        updated_mutafy_refseq_warnings.sort_values(by='gene', inplace=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_refseq_warnings.to_csv(f'{mutafy_directory}/05_refseq_warnings_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the data for all the structures for all genes of the current webrun
        refseq_warnings = updated_mutafy_refseq_warnings[updated_mutafy_refseq_warnings.gene.isin(list(fasta_df.gene_name.unique()))]
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        refseq_warnings.to_csv(f'{mutafy_directory}/05_refseq_warnings_mutafy.csv', index = False)
    
    if exists(f'{mutafy_directory}/05_blastp_warnings_mutafy.csv'):
        mutafy_blastp_warnings = pd.read_csv(f'{mutafy_directory}/05_blastp_warnings_mutafy.csv')
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        # in order to do that, we convert all the values in the refseq_warnings df to strings
        blastp_warnings = blastp_warnings.astype('str')
        # now we can concatenate the two dfs
        # we update the mutafy data, concatenate it with our df and drop potential duplicates
        updated_mutafy_blastp_warnings = pd.concat([mutafy_blastp_warnings, blastp_warnings], ignore_index=True).drop_duplicates()
        # we sort the df again first according to gene name 
        updated_mutafy_blastp_warnings.sort_values(by='gene', inplace=True)
        # we write the updated mutafy data to a csv file (which overwrites the one from the previous mutafy run)
        updated_mutafy_blastp_warnings.to_csv(f'{mutafy_directory}/05_blastp_warnings_mutafy.csv', index=False)        
        # we get a slice of the mutafy data with all the genes of the current search and save this to the results folder
        # this is a df with all the data for all the structures for all genes of the current webrun
        blastp_warnings = updated_mutafy_blastp_warnings[updated_mutafy_blastp_warnings.gene.isin(list(fasta_df.gene_name.unique()))]
    else:
        # if the mutafy file doesn't exist yet, we write it out with the data from the current run
        blastp_warnings.to_csv(f'{mutafy_directory}/05_blastp_warnings_mutafy.csv', index = False)

# we change back to the results directory
os.chdir(results_dir)

# write csv output files
print(f'>>> writing csv file containing blastp results for all structures of all genes\n')
xml_df.to_csv('05_blastp_results.csv', index = False)
# we also write the warnings to csv files
refseq_warnings.to_csv('05_refseq_warnings.csv', index = False)
blastp_warnings.to_csv('05_blastp_warnings.csv', index = False)
                
# Now if this is a webrun, we will also update the first mutafy csv file (00_search_overview_PDBids_mutafy.csv)
# as all the new pdb ids associated with the current run have now been downloaded/parsed and blasted,
# so we don't need to do this again in a future run!
# (Potential Problem: except if we want to change some of the parameters in the earlier scripts, e.g. blast parameters, but
# we don't do this for now, i.e. currently no way to distinguish between different cutoffs used in previous runs..)
# we read in the mutafy data (file called 00_search_overview_PDBids_mutafy.csv) from the mutafy_directory
# and the file of the current run called 00_search_overview_PDBids.csv from the results_dir
# then we will combine the two / updated the mutafy data
# we have to use the rows from the data for the current run and overwrite the corresponding rows
# in the mutafy data to do this
# read in the data
webrun_data = pd.read_csv(f'{results_dir}/00_search_overview_PDBids.csv', usecols=['gene_name', 'n_available_structures', 'available_structures'])

# for the very first webrun (or whenever we delete this file for some reason), the mutafy file won't exist
# so we add a try and exept statement and if the file doesn't exist, we write the webrun data to a csv file as a first
# mutafy data file
try:
    mutafy_data = pd.read_csv(f'{mutafy_directory}/00_search_overview_PDBids_mutafy.csv', usecols=['gene_name', 'n_available_structures', 'available_structures'])
except FileNotFoundError:
    mutafy_data = None
    webrun_data.to_csv(f'{mutafy_directory}/00_search_overview_PDBids_mutafy.csv', index=False)
    
# now we update the mutafy_data df if it exists
# we can update it with the webrun_data using DataFrame.update, which aligns on indices
# (https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.update.html)
if mutafy_data is not None:    
    mutafy_data.set_index('gene_name', inplace=True) # we need the gene name to be the index
    mutafy_data.update(webrun_data.set_index('gene_name')) # we update the df with data from the webrun    
    mutafy_data.sort_values(by='gene_name', inplace=True, ignore_index=True) # we sort the df 
    mutafy_data.reset_index(inplace=True)  # reset index
    # now that we updated the mutafy data, we write it to a csv file to be used for future runs
    # all structures listed in this file won't be downloaded/parsed/blasted again when running mutafy :-)
    mutafy_data.to_csv(f'{mutafy_directory}/00_search_overview_PDBids_mutafy.csv', index=False)  
        
# change back to target directory
os.chdir(target_directory)
    
print('\n============================== Summary ================================================\n')
if web_run:
    print(f'Complete! \n    Performed BLASTp on a total of {len(xml_slice_to_add)} sequences for {len(xml_slice_to_add.gene_name.unique())} gene(s) listed in the file 04_fasta_combined_info.csv.\n')

elif not web_run:
    print(f'Complete! \n    Performed BLASTp on a total of {len(fasta_df)} sequences listed in the file 04_fasta_combined_info.csv.\n')
    print('     All reference sequences (one per gene) used for blastp have been stored in .fasta format in the Results/RefSeqs folder.')
    print('     All all unique sequences per structure used for blastp have been stored in .fasta format in the Results/RefSeqs/PDB_seqs_and_blastp_outputs folder.')
    print('     All BLASTp outputs are stored in .xml format in the Results/RefSeqs/PDB_seqs_and_blastp_outputs folder.')

print('\nThe following files have been created and stored in the Results folder:')
print('   o      05_blastp_results.csv                (lists all info contained in the input file as well as their blastp results)')
print('   o      05_refseq_warnings.csv                (lists genes with no or more than one identified reference sequence (only first one is used for further analyses))')
print('   o      05_blastp_warnings.csv                (lists warnings regarding blastp, including if blastp failed or if there is more than one alignment (should only be one as only one reference is used)\n\n')

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