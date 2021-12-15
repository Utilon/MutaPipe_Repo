# This script takes a csv file (03_fasta_combined_info.csv) containing information for all genes extracted from their
#  respective fasta files and fasta_ex files as input and will:
#      - download the reference sequence for the gene from uniprot
#      - perform BLASTp with each sequence of each structure (=each row of the df) against the respective
#        reference sequence (e.g. FUS canonical sequence for all sequences in all FUS structures)
#      - uses the blast output (xml files) to identify mismatches and add this information to the df
#      - outputs the following files:
#                - all reference sequences and unique sequences of all structures are stored in .fasta format in the newly created Results/RefSeqs folder
#                - all BLASTp outputs are stored in .xml format in the newly created Results/RefSeqs folder
#                - a csv file called 04_blast_two_sequences.csv listing all the information in the input file and the corresponding blastp results

# ===========================================================================================

# Set up
import pandas as pd
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import re
import sys
import argparse
from datetime import datetime

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

uniprot_fasta = f'{target_directory}/Uniprot_reference_seqs/UP000005640_9606.fasta' # specify path to uniprot reference fasta
                                            
                                            
# Now we create an argument parser called ap to which we can add the arguments we want to have in the terminal
ap = argparse.ArgumentParser(description="""****    This script takes a csv file (03_fasta_combined_info.csv) containing information for all input genes extracted from their
#  respective fasta files and fasta_ex files as input and will: 
1. retrieve the reference sequence for the gene associated with each structure from the uniprot reference fasta file 
2. perform BLASTp with each sequence of each structure against the respective reference sequence (e.g. FUS canonical sequence for all sequences in all FUS structures) 
3. uses the blast output (xml files) to identify mismatches / variants in all sequences 
4. outputs the following files: 
(1) all reference sequences and unique sequences of all structures are stored in .fasta format in the newly created Results/RefSeqs folder 
(2) all BLASTp outputs are stored in .xml format in the newly created Results/RefSeqs folder 
(3) a csv file called 04_blast_two_sequences.csv listing all the information on all sequences in all structures for all genes and the corresponding blastp results    ***""")

ap.add_argument('-bp','--blastp_path', required=False, help=f"Specify the path to blastp on your system ; default = {blastp_path}")
ap.add_argument("-l", "--log", type=str2bool, required = False, help=f'Write output to .log file in output directory if set to True, default = {str(create_search_log)}')
ap.add_argument("-t", "--target", required = False, help=f'Specify target directory, default = {target_directory}')
ap.add_argument("-refseq", "--reference_sequences", required = False, help=f'Specify path to uniprot reference fasta, default = {uniprot_fasta}')

args = vars(ap.parse_args())

# Now, in case an argument is used via the terminal, this input has to overwrite the default option we set above
# So we update our variables whenever there is a user input via the terminal:
blastp_path = blastp_path if args["blastp_path"] == None else args["blastp_path"]
create_search_log  = create_search_log  if args["log"]   == None else args["log"]
target_directory  = target_directory if args["target"]   == None else args["target"]
uniprot_fasta = uniprot_fasta if args["reference_sequences"] == None else args["reference_sequences"]

# ----------------------------------------------------------------------------------------------------------------------------------

# We want to write all our Output into the Results directory

results_dir = f'{target_directory}/Results' #define path to results directory

# ----------------------------------------------------------------------------------------------------------------------------------

#  create log file for console output:
if create_search_log == True:
    with open(f'{results_dir}/search_log_04.txt', 'w') as search_log:
        search_log.write(f'Search log for 04_blast_against_reference.py\n\n')
    sys.stdout = open(f'{results_dir}/search_log_04.txt', 'a')

# store current date and time in an object and print to console / write to log file
start_time = datetime.now()
print(f'start: {start_time}\n')


# ----------------------------------------------------------------------------------------------------------------------------------


# we make a new folder in the Results folder to store all reference sequences (permanently) and individual fasta files (temporarily)
# which we need during our for loop when doing the blastp
refseqs_dir = f'{results_dir}/RefSeqs'
if not os.path.isdir(refseqs_dir):
    os.mkdir(f'{results_dir}/RefSeqs')
    
# change to refseqs_dir
os.chdir(refseqs_dir)

# Read in data from csv files
fasta_df = pd.read_csv(f'{results_dir}/03_fasta_combined_info.csv')


# create a df based on the fasta_df with additional columns to populate with results from blast p
xml_df = pd.concat([fasta_df, pd.DataFrame(columns=['alignment_length', 'hsp_number',
                                                    'hsp_length', 'score', 'bit', 'e-value', 'similarity',
                                                    'mismatches_incl_gaps', 'mismatches_excl_gaps',
                                                    'close_mismatches', 'mismatch_substitutions',
                                                    'close_mismatch_substitutions',
                                                    'gaps', 'gaps_in_query', 'gaps_in_sbct', 'gaps_start_pos_query',
                                                    'gaps_length_query', 'query_del_sbcjt_ins', 'gaps_start_pos_sbjct',
                                                    'gaps_length_sbjct', 'sbjct_del_query_ins', 'query_sequence',
                                                    'sbjct_sequence', 'match_sequence'])])                                            

# we also create empty dfs to capture all our WARNINGS (as defined by me, haha) and write them to a csv file at the end
# a df to store all genes for which no or more than one reference sequence is found in the uniprot reference fasta
refseq_warnings = pd.DataFrame(columns=['gene', 'n_refseqs_identified'])
blastp_warnings = pd.DataFrame(columns=['gene', 'blast_status', 'n_alignments'])

# we also read in the a fasta file downloaded from uniprot containing all canonical sequences for the human genome:
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
        print(f'>>> Writing canonical reference sequence for {gene} (gene {gene_counter} of {len(fasta_df.gene_name.unique())}) to fasta file for blastp')
        # we make a list and include only sequences/fasta records that contain the gene name in the description
        # in order to get the correct sequence, e.g. NEK1 instead of NEK10, we include a space after the gene name!!!
        potential_reference_sequences = [sequence for sequence in uniprot_refseqs if f'GN={gene + " "}' in sequence.description]
        # there should only be one reference sequence, so we take the first and only element of this list as our reference sequence
        # but e.g. for SMN1 there's no reference at all (its name is snm2 for some reason!? download new uniprot fasta!)
        # so we add a try and and except statement:
        try:
            reference_sequence = potential_reference_sequences[0]
        except IndexError:
            print(f'WARNING!        No reference sequence found for gene {gene}')
            refseq_warnings.loc[len(refseq_warnings)] = [gene, 0]
            continue
        # update previous_refseq with the current reference sequence
        previous_refseq = reference_sequence
        # we write a fasta file for the current reference sequence that we can later use to blastp
        # we define a variable for the name of this file, so we can read it in again later
        refseq_fasta = f'{gene}_reference.fasta'
        # now we write the fasta file 
        SeqIO.write(reference_sequence, refseq_fasta, 'fasta')
        # we update the variable previous_refseq_fasta and the previous_gene variable:
        previous_gene = gene
        previous_refseq_fasta = refseq_fasta
        # print WARNING if more than one reference sequence has been identified (in this case the first one is taken as a reference sequence for the next step)
        if len(potential_reference_sequences) > 1:
            print(f'WARNING: more than 1 reference sequence identified for gene {gene}: {len(potential_reference_sequences)} potential reference sequences')
            refseq_warnings.loc[len(refseq_warnings)] = [gene, len(potential_reference_sequences)]
    # if this gene is the same as the previous on, the current reference sequence is the same as the previous one too
    elif gene == previous_gene:
        # the reference sequence is the same as the previous one, so we don't need to write a new fasta file (it already exists)
        # in order to prevent potential bugs, we still update our variable reference_sequence for this gene,
        # even if this is potentially redundant as we don't need this variable but the fasta file
        reference_sequence = previous_refseq
        # the fasta file is the same as the previous one
        refseq_fasta = previous_refseq_fasta
        
    # so now that we have the sequence and the reference sequence in fasta format,
    # we want to perform a blastp of the two
    # useful info on how to do this here: https://www.biostars.org/p/42687/
    # also here: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec127
    
    # Run BLAST and save output as xml file
    print(f'        >>> perfoming blastp for {gene} fasta file {fasta}')
    blastp_cline = NcbiblastpCommandline(cmd=blastp_path, query=fasta, subject= refseq_fasta, outfmt=5, out=fasta.replace('fasta', 'xml'))
    try:
        stdout, stderr = blastp_cline()
    except:
        print(f'        blastp UNSUCCESFUL for {gene} fasta {fasta}')
        blastp_warnings.loc[len(blastp_warnings)] = [gene, f'FAILED for {fasta}', np.nan]
        continue    
    
    # parse xml file/blastp output to read blastp results
    print(f'        >>> parsing blastp output for fasta file {fasta}\n')
    result_handle = open (fasta.replace('fasta', 'xml'))
    blast_record = NCBIXML.read(result_handle)
    
    # now we can use the blast_record to extract information on mismatches, gaps etc.
    #there is only one alignment in the blast_record (only one subject sequence given, i.e. our reference sequence for each gene)
    # or there could be no alignments, so we add a try and except statement
    try:
        alignment = blast_record.alignments[0]
    except:
        blastp_warnings.loc[len(blastp_warnings)] = [gene, f'No alignments for {fasta}', 0]
        continue
    
    # in case there are more than one alignment: print warning
    if len(blast_record.alignments) > 1:
        print(f'WARNING: more than one alignment found in blastp output {fasta.replace("fasta", "xml")}')
        print('            NOTE: only first top alignment used in further analysis')
        blastp_warnings.loc[len(blastp_warnings)] = [gene, f'>1 alignment for {fasta}', len(blast_record.alignments)]
        
    hsp_counter = 0
    # there may be multiple high scoring segment pairs (HSPs) in each alignment, if that's the case, we only take the first one for now
    # (otherwise, we'd have to create a new row for each hsp in the output df or add more columns for more hsp data to the df... currently,
    # if we delete the [:1] and loop over all hsps instead, only the data of the last hsp get appended to the output df)
    # for hsp in alignment.hsps[:1]:
    # better idea: we choose the best hsp (with the highest bit score) and save it in a variable called hsp
    hsp = None
    best_bit_score = 0
    best_hsp_counter = 0 
    for this_hsp in alignment.hsps:
        hsp_counter += 1
        if this_hsp.bits > best_bit_score:
            best_bit_score = this_hsp.bits
            hsp = this_hsp
            best_hsp_counter = hsp_counter
                  
    # find the number of mismatches and close mismatches using the hsp.match sequence
    mismatches = len([letter for letter in hsp.match if letter == ' '])
    close_mismatches = len([letter for letter in hsp.match if letter == '+'])

    # find the positions of mismatches and close mismatches in the hsp.match sequence
    mismatch_pos = [char for char, ltr in enumerate(hsp.match) if ltr == ' ']
    close_mismatch_pos = [char for char, ltr in enumerate(hsp.match) if ltr == '+']

    # define two empty lists to populate with all mismatch and close mismatch AA substitutions                                                   
    mismatch_substitutions = []
    close_mismatch_substitutions = []
    
     # =======================================================================
            
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
    # create empty dict to populate with {'query_seq' : 'sbjct_seq', ...}
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
    # create empty dict to populate with {'sbjct_seq' : 'query_seq', ...}
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
    
    # now we have all the information from the blastp output stored in variables and can store them in a dataframe:
    xml_df.loc[fasta_df_index, 'alignment_length'] = alignment.length
    xml_df.loc[fasta_df_index, 'hsp_number'] =  best_hsp_counter     
    xml_df.loc[fasta_df_index, 'hsp_length'] = hsp.align_length
    xml_df.loc[fasta_df_index, 'score'] = hsp.score
    xml_df.loc[fasta_df_index, 'bit'] = hsp.bits
    xml_df.loc[fasta_df_index, 'e-value'] = hsp.expect
    xml_df.loc[fasta_df_index, 'similarity'] = hsp.identities / hsp.align_length
    xml_df.loc[fasta_df_index, 'mismatches_incl_gaps'] = mismatches
    xml_df.loc[fasta_df_index, 'mismatches_excl_gaps'] = mismatches_excl_gaps
    xml_df.loc[fasta_df_index, 'close_mismatches'] = close_mismatches
    xml_df.loc[fasta_df_index, 'mismatch_substitutions'] = mismatch_substitutions
    xml_df.loc[fasta_df_index, 'close_mismatch_substitutions'] = close_mismatch_substitutions
    xml_df.loc[fasta_df_index, 'gaps'] = hsp.gaps
    xml_df.loc[fasta_df_index, 'gaps_in_query'] = all_gaps_query
    xml_df.loc[fasta_df_index, 'gaps_in_sbct'] = all_gaps_sbjct
    xml_df.loc[fasta_df_index, 'gaps_start_pos_query'] = gaps_start_pos_query
    xml_df.loc[fasta_df_index, 'gaps_length_query'] = gaps_len_query
    xml_df.loc[fasta_df_index, 'query_del_sbcjt_ins'] = str(gap_seqs_dict_query)
    xml_df.loc[fasta_df_index, 'gaps_start_pos_sbjct'] = gaps_start_pos_sbjct
    xml_df.loc[fasta_df_index, 'gaps_length_sbjct'] = gaps_len_sbjct
    xml_df.loc[fasta_df_index, 'sbjct_del_query_ins'] = str(gap_seqs_dict_sbjct)
    xml_df.loc[fasta_df_index, 'query_sequence'] = hsp.query
    xml_df.loc[fasta_df_index, 'sbjct_sequence'] = hsp.sbjct
    xml_df.loc[fasta_df_index, 'match_sequence'] = hsp.match


# Before we change back to the results folder, we clean up the RefSeqs folder a bit
# we want to make a subdirectory in this folder called PDB_seqs_and_blastp_outputs
# we will move all the fasta files from the PDB structures and the blastp outputs there and
# only keep the reference fasta directly in the RefSeqs directory
# first we get a list of all files in the current folder (RefSeqs)
all_fasta_files_for_blastp = [f for f in listdir(refseqs_dir) if isfile(join(refseqs_dir, f))]

# now we get a list of all the files that are not reference fasta files and move them to a new subfolder
PDB_fasta_files_for_blastp = [f for f in all_fasta_files_for_blastp if not 'reference.fasta' in f]
# make the folder
os.mkdir('PDB_seqs_and_blastp_outputs')
# and move the files
for PDB_fasta in PDB_fasta_files_for_blastp:
    os.rename(f'{PDB_fasta}', f'PDB_seqs_and_blastp_outputs/{PDB_fasta}')
              
# we change back to the results directory
os.chdir(results_dir)

# write csv file
print(f'>>> writing csv file containing blastp results for all structures of all genes\n')
xml_df.to_csv('04_blast_two_sequences.csv', index = False)

# we also write the warnings to csv files
refseq_warnings.to_csv('04_refeseq_warnings.csv', index = False)
blastp_warnings.to_csv('04_blastp_warnings.csv', index = False)
    
# change back to target directory
os.chdir(target_directory)
    
    
print('\n============================== Summary ================================================\n')
print(f'Complete! \n    Performed BLASTp on a total of {len(fasta_df)} sequences listed in the file 03_fasta_combined_info.csv.\n')

print('     All reference sequences (one per gene) used for blastp have been stored in .fasta format in the Results/RefSeqs folder.')
print('     All all unique sequences per structure used for blastp have been stored in .fasta format in the Results/RefSeqs folder.')
print('     All BLASTp outputs are stored in .xml format in the Results/RefSeqs folder.')

print('\nThe following files have been created and stored in the Results folder:')
print('   o      04_blast_two_sequences.csv                (lists all info contained in the input file as well as their blastp results)')
print('   o      04_refeseq_warnings.csv                (lists genes with no or more than one identified reference sequence (only first one is used for further analyses))')
print('   o      04_blastp_warnings.csv                (lists warnings regarding blastp, including if blastp failed and if there are more than one alignment (should only be one as only one reference is used)\n\n')

# store current date and time in an object and print to console / write to log file
end_time = datetime.now()
print(f'start: {start_time}\n')
print(f'end: {end_time}\n\n')

# close search log
if create_search_log == True:
    sys.stdout.close()
        