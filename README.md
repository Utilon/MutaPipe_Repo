# :dna::dna::dna:    **M    u    t    a    P    i    p    e**     :dna::dna::dna:   
## A bioinformatics pipeline to identify high quality mutant and WT PDB structures for genes of interest
```diff
+TO BE NOTED: We are always working to improve MutaPipe, so any bug reports or suggestions are highly welcome.
```
## Table of Contents
1. [Introduction](#introduction)
2. [Citation](#citation)
3. [Documentation](#documentation)
  	* [Workflow](#workflow)
  	* [Incorporated scripts](#incorporated-scripts)
	* [Minimum Requirements](#minimum-requirements)
	* [Local Deployment](#local-deployment)
	* [Usage](#usage)
	* [Usage Example](#usage-example)
	* [Output](#output)
	* [How to Download the Reference Proteome](#how-to-download-the-reference-proteome)

## Introduction

MutaPipe is a fast and efficient bioinformatics pipeline to screen the Protein Data Bank (PDB) for genes of interest and retrieve the highest quality protein structure (best resolution) for each unique sequence associated with a gene (currently only available for the species Homo sapiens). Additionally, whenever corresponding data is available, variants will be annotated using information on variant pathogenicity from ClinVar.

This allows researchers to efficiently screen the PDB for the most suitable template structures for a specific WT or mutant gene/protein which can be used for further in silico analysis like mutagenesis experiments or molecular dynmamics simulations.

![alt text](https://github.com/Utilon/Pipeline_Git/blob/main/files_for_README/Figure1_pipeline_simple_flowchart.jpeg)

Figure 1. *Pipeline overview. MutaPipe requires a list of gene names as input. The pipeline firstly performs a search of the entire Protein Data Bank (PDB) to retrieve information on all protein structures available in the PDB which are associated with the input genes (and the species: Homo sapiens). In a second step, a BLASTp against the canonical reference sequence for a given input gene (obtained from Uniprot) is performed for each sequence of every identified protein structure associated with said gene of interest. The BLASTp alignment is used to identify mismatches in the amino acid sequence of the PDB structures relative to the canonical sequence. In a subsequent step, the data is filtered and rearranged to output a table listing the highest quality structure for all available sequences in the PDB which are associated with the input genes. Additionally, if a given amino acid change in an input gene is listed in ClinVar, this information will be added to the output table.*



## Citation

MutaPipe has been developed by Deborah Ness, an NIHR Maudsley PhD student, in collaboration with and supervised by Dr Alfredo Iacoangeli, for her PhD project at the Department of Basic and Clinical Neuroscience at King's College London.

This work has not yet been published (Dec 2021). 

However, a dataset produced with Mutapipe for all genes available in the ClinVar Database has been presented as a poster at an international research conference in Italy in September 2021 ("THE 1° MASBIC - DISVA ANNUAL SYMPOSIUM - From structure to function: unveiling the role of proteins in health and disease")
The dataset is available [here](https://drive.google.com/drive/folders/1hqvAOLGMJuYqKuB1poH5_JPXFHDxwX80?usp=sharing).

## Documentation

### Workflow 

MutaPipe uses the inputted genes names (also possible in .txt format) to find all associated PDB IDs (using the PDB API) and download all corresponding fasta, mmCIF, and PDB files for further anaylses (will need storage depending on how many structures are available for the input genes (and their size)).
In a next step, MutaPipe parses all the downloaded files and combines all relevant information (e.g. resolution, sequence, unsolved residues in PDB structures). Subsequently, MutaPipe accessess the canonical protein sequence for all input genes (from Uniprot fasta file) and uses BLASTp (e-value threshold: .05) to identify all variants available in the identified PDB structures. BLASTp is performed for every sequence in every structure against the Uniport reference sequence for the gene associated with this structure (e.g. BLASTp of all sequences in a given FUS structure against the canonical sequence for FUS). The .xml output files generated by BLASTp are parsed and used to extract information, incl.: e-value, bit, alignment length, mismatches, close mismatches, gaps, indels.  This information will be used by MutaPipe identify all unique sequences in the PDB which are associated with the input genes and subsequently get the best structure (the one with the highest resolution) for all for each of them. 

MutaPipe achieves this by sorting the data to get:
- the best structure per point mutation (structures with only this one mutation and no other mutations)
- the best structure per unique combination of mutations available in the PDB (including WT/no mutation)
- the best structure for any specific mutation, regardless of other mutations in the same structure

Finally, whenever corresponding data is available, MutaPipe will annotate variants using variant information from ClinVar, incl. information on variant pathogenicity.

**Filtering**
To ensure sequences not associated with the gene of interest which are present in a PDB structure (e.g. if a protein of interest is co-crystallised with another peptide) do not get included in the final output, MutaPipe also has an additional step to filter out sequences:
1. which are shorter than a given percentage of the reference sequence (set variable relative_sequence_length; default value = 0.5)
2. whose best hsp covers less than a given percentage of the reference sequence (set variable hsp_coverage; default value = 0.1)




![alt text](https://github.com/Utilon/Pipeline_Git/blob/main/files_for_README/Figure2_Pipeline_detailled_workflow.jpeg)

Figure 2. *Detailled description of the steps performed by MutaPipe to identify variant and wildtype PDB structures associated with the input genes and identify the highest quality structure for each available sequence. (Note: The final step which adds available variant annotations from ClinVar is not displayed).*


### Incorporated scripts

There are currently 9 different python scripts which are incorporated in MutaPipe:

| Script                                                                                                                        | Input                                                                                                                                                                                   | Operations                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         | Output                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|-------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [00](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/00_search_pdb.py)                                              | genes in .txt format<br>(e.g. SOD1 ALS2 FUS)                                                                                                                                            | - creates a folder called Results in the current <br>working directory where all the output is going to be stored<br><br>- searches the pdb for all structures associated with <br>each gene name (in Homo Sapiens)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | - directory: Results/<br><br>- `00_search_overview_PDBids.csv`<br>*contains all gene names and corresponding PDB IDs if available*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| [01](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/01_download_files.py)                                          | `00_search_overview_PDBids.csv`                                                                                                                                                         | - creates a folder for each gene in the Results directory<br><br>- downloads mmCIF, pdb, and fasta files for all structures into the respective gene folder                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | - directory per gene: Results/GENENAME/<br><br>In each respective gene folder:<br>- mmCIF, pdb, and fasta files<br>stored in format `pdbID.cif/pdb/fasta`<br><br>In the Results folder:<br>- `01_search_overview_folders.csv`<br>*lists all the the newly created folders and their contents*<br>- `01_search_overview_n_structures.csv`<br>*lists number of structures retrieved per gene*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| [02](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/02_parse_cif_files.py)                                         | `01_search_overview_folders.csv`                                                                                                                                                        | **loops over all created folders and parses all mmCIf files, extracting:**<br><br>- resolution (999 for missing values/NMR structures)<br>- polypeptide sequences (corresponds to sequences as shown in PyMOL)<br>- fasta sequences                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | In each respective gene folder:<br>- `_ex.fasta` file for every structure <br>*fasta file created from the mmCIF file*<br>- `GENENAME_02_resolutions.csv`<br>*contains the resolution for all structures associated with this gene*<br>- `GENENAME_02_poly_seq.csv`<br>*contains all polypeptide sequences for all structures associated with this gene*<br><br>In the Results folder:<br>- `02_all_resolutions.csv`<br>*contains the resolutions of all parsed structures for all genes*<br>- `02_all_poly_seq.csv`<br>*contains all polypeptide sequences for all structures of all genes*                                                                                                                                                                                                                                                                                                                                                                 |
| [03](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/03_parse_fasta_files.py)                                       | `01_search_overview_folders.csv`                                                                                                                                                        | **loops over folders containing fasta and _ex.fasta files and:**<br><br>extracts info from fasta files, incl:<br>   - chain name<br>   - description<br>   - species<br>   - sequence<br><br>extracts info from _ex.fasta files, incl:<br>   - chain name<br>   - description<br>   - uniprot id <br>   - sequence<br><br>combines info from the two files                                                                                                                                                                                                                                                                                                                                                                                                                                                         | In each respective gene folder:<br>- `GENENAME_03_fasta_info.csv`<br>*contains information extracted from all fasta files for this gene*<br>- `GENENAME_03_fasta_ex_info.csv`<br>*contains information extracted from all _ex.fasta files for this gene*<br>- `GENENAME_03_fasta_combined_info.csv`<br>*contains combined information extracted from all fasta and _ex.fasta files for this gene*<br><br>In the Results folder: <br>- `03_fasta_info.csv`<br>*contains information extracted from all fasta files for all genes*<br>- `03_fasta_ex_info.csv`<br>*contains information extracted from all _ex.fasta files for all genes*<br>- `03_fasta_combined_info.csv`<br>*contains combined information extracted from all fasta and _ex.fasta files for all genes*                                                                                                                                                                                      |
| [04](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/04_blast_against_reference.py)                                 | `03_fasta_combined_info.csv`                                                                                                                                                            | - creates a directory called RefSeqs in the Results directory<br>............................................................................................................<br><br>- loops over input csv file (one row for each unique sequence in all the pdb structures for all genes)<br><br>- extracts the reference sequence for each gene from the Uniprot reference fasta into the RefSeqs directory<br><br>- writes a fasta file for each unique sequence/chain in the input csv into the RefSeqs directory<br><br>- performs BLASTp of all sequences against their respective reference sequence (e.g. FUS canonical sequence serves as reference for all sequences in all FUS structures) (output stored in .xml format in RefSeqs)<br><br>- uses the blast output (xml files) to identify mismatches | - directory: Results/Refseqs<br><br>In the RefSeqs folder:<br>- reference sequence fasta files<br>in format `GENE_reference.fasta`<br>- unique sequences fasta files<br>in format `GENE_pdbID_Chains.fasta`<br>- BLASTp output files (.xml)<br>in format `GENE_pdbID_Chains.xml`<br><br>In the Results folder: <br>- `04_blast_two_sequences.csv`<br>*lists all the information in the input file and the corresponding BLASTp results*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| [05](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/05_pdb_extract_unsolved_res.py)                                | `01_search_overview_folders.csv`<br><br>`04_blast_two_sequences.csv`                                                                                                                    | - loops over folders and extracts info on missing residues / residues which have not been solved in the crystal structure from each pdb file<br><br>- combines information on unsolved residues with info from 06_blast_fasta.csv / 08_blast_two_sequences.csv                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | In the Results folder: <br>- `05_unsolved_residues_per_structure.csv`<br>*lists all unsolved residues in all structures for all genes (one row for each structure)*<br>- `05_unsolved_residues_per_chain.csv`<br>*lists all unsolved residues in all chains of all structures for all genes (one row for each chain)*<br>- `05_all_info.csv`<br>*contains all the information from 04_blast_two_sequences.csv plus on unsolved residues extracted from pdb files*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| [06](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/06_best_structure_per_mutation.py)                             | `05_all_info.csv`<br><br>`02_all_resolutions.csv`                                                                                                                                       | - combines the two dfs (according to PDBid)<br><br>- filters out sequences shorter than a given percentage of the reference sequence (set variable `relative_sequence_length`)<br><br>- filters out sequences whose best hsp covers less than a given percentage of the reference sequence (set variable `hsp_coverage`)<br><br>- sort/filter the df to get:<br> 1. the best structure per point mutation (structures with only this one mutation and no other mutations)<br> 2. the best structure per unique combination of mutations available in the PDB<br> 3. the best structure for any specific mutation, regardless of other mutations in the same structure                                                                                                                                              | In each respective gene folder:<br>- `GENENAME_06_best_structure_per_point_mutation.csv`<br>*lists best structure for each point mutation (one mutation per structure) in this gene*<br>- `GENENAME_06_best_structure_all_unique_combinations.csv`<br>*lists best structure for all unique mismatch combinations for this gene*<br>- `GENENAME_06_best_structure_any_mutation.csv`<br>*lists best structure for any mismatch in this gene regardless of other mismatches in this structure*<br><br>In the Results folder: <br>- `06_best_structure_per_point_mutation.csv`<br>*lists best structure for each point mutation (one mutation per structure) in all genes*<br>- `06_best_structure_all_unique_combinations.csv`<br>*lists best structure for all unique mismatch combinations for all genes*<br>- `06_best_structure_any_mutation.csv`<br>*lists best structure for any mismatch for all genes regardless of other mismatches in this structure* |
| [07_a](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/07_a_ClinVar_Annotations_edirect_per_gene_download_files.py) | `00_search_overview_availability.csv`                                                                                                                                                   | - query ClinVar for information on each gene of interest and download xml files with ids for all variants using edirect<br><br>- parse xml files with ClinVar ids and download/retrieve xml data from ClinVar for all variant ids associated with a gene                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | - all xml files downloaded from ClinVar are stored in the newly created ClinVar_Annotations folder (a subfolder in the Results directory)<br><br>- `07_a_ClinVar_Annotations_genes_no_data_retrieved.txt`<br>*lists all genes for which no ClinVar annotations could be retrieved*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| [07_b](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/07_b_ClinVar_Annotations_edirect_per_gene_parse_files.py)    | `00_search_overview_availability.csv`                                                                                                                                                   | - parse xml files and create a df with ClinVar information for all variants for all input genes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | - `07_b_ClinVar_Annotations.csv`<br>*lists ClinVar annotations for all variants in all input genes*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| [08](https://github.com/Utilon/Pipeline_Git/blob/main/MutaPipe/08_add_clinvar_annotations_to_best_structures.py)             | `06_best_structure_all_unique_combinations.csv`<br><br>`06_best_structure_any_mutation.csv`<br><br><br>`06_best_structure_per_point_mutation.csv`<br><br>`07_b_ClinVar_Annotations.csv` | - adds availalable ClinVar annotations to all three best_structure tables (csv-files)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | - `08_best_structure_all_unique_combinations.csv`<br>*lists best structure for all unique mismatch combinations for all genes (incl. ClinVar annotations)*<br>- `08_best_structure_any_mutation.csv)`<br>*lists best structure for any mismatch for all genes regardless of other mismatches in this structure(incl. ClinVar annotations*<br>- `08_best_structure_per_point_mutation.csv`<br>*lists best structure for each point mutation (one mutation per structure) in all genes (incl. ClinVar annotations)*                                                                                                                                                                                                                                                                                                                                                                                                                                            |


### Minimum Requirements (incomplete)

- python > 3.6.8 
- RAM: ??
- Space required by the installation: 36 MB
- Scratch space for usage: depends the number and and size of PDB structures associated with the input genes. To check how many structures are available for your genes of interest, run `python3 00_search_pdb.py -g GENES_OF_INTEREST`. This will generate console output as well as csv files both of which provide an overview of how many structures have been identified for your input genes.

### Local Deployment

To obtain MutaPipe please use git to download the most recent development tree:

```bash
git clone https://github.com/Utilon/Pipeline_Git.git
```
#### NOTE: Reference Proteome Required
In order to run MutaPipe a fasta file containing the reference proteome / canonical protein sequences for the organism in question (= Homo sapiens) will have to be downloaded from Uniprot and stored in a folder called `Uniprot_reference_seqs` in the `MutaPipe` directory (where the python scripts are stored). The filename should be `UP000005640_9606.fasta` for the script to work properly.
The latest version of this fasta file can be accessed via the [Uniprot ftp client](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/).

**NOTE**: When installing MutaPipe a fasta file of the reference proteome for Homo sapiens will be included (version from December 2021). Please feel free to update to a newer version if one is available (we will also regularly update the reference proteome file on here, but this is currently not done automatically).


### Usage (incomplete)

Its basic use requires the following options:

```bash

  -format FORMAT        options are bam, sam, fastq, vcf [string] 
  -out OUT              path to the output folder. It has to end in /" e.g. /home/user/local/test_folder/

 ```
 The desired pipeline stages are performed according to the optional arguments selected:
 
 ```bash
  -filter_string FILTER_STRING  bcftools filter string, eg GQ>20 & DP>10 (Default = "")
  -iobio                if this flag is set iobio services will be started at the end of the analysis (Default = "False")
```

Finally, a set of optional arguments can be used to customise the analysis:

 ```bash
-RG RG                if this flag is set the alignment stage will add the provided in paths_and_configs.py read group (Default = "False")
-paired PAIRED        options are 1 for paired end reads and 0 for single end reads (Default = "1")
```

#### Usage Example (incomplete)

Let's assume we have human paired end whole exome sequening data in two fastq files and want to perform snvs/indels calling vs hg19, annotation and explore the results using the iobio services. The DNAscan command line would be:

 ```bash
python3 /path/to/DNAscan/scripts/DNAscan.py -format fastq -in data1.fq.gz -in2 data2.fq.gz -reference hg19 -alignment -variantcalling -annotation -iobio -out /path/to/outdir/ -mode fast
```
Using the sequencing data provided in the data folder:

 ```bash
cd /path/to/DNAscan_main_dir
 
python3 scripts/DNAscan.py -format fastq -in data/test_data.1.fq.gz -in2 data/test_data.2.fq.gz -reference hg19 -alignment -variantcalling -annotation -iobio -out outdir/ -mode fast -BED
```
IMPORTANT: All paths in DNAscan end with "/"


### Output (incomplete)

DNAscan output tree:

```bash

./$out_dir-| # this is the folder given to DNAsca using the -out flag. It will contain the aligned sequecing data ($sample_name.bam) as well as some temporanery files
           |
           |-results # This will contain the output of the analyses. E.g. $sample_name_sorted.vcf.gz , $sample_name_SV.vcf.gz, virus_results.txt, etc  
           |
           |-reports # If any report flags is used, this folder will contain the reports. E.g. $sample_name_vcfstats.txt if the -calls_report flag is used
           |
           |-logs # Logs files are generated in this folder
 
```          
           

### How to Download the Reference Proteome

MutaPipe will need a fasta file with all canonical protein sequences of the species of interest.

In order to run MutaPipe a fasta file containing the reference proteome / canonical protein sequences for the organism in question (= Homo sapiens) will have to be downloaded from Uniprot and stored in a folder called `Uniprot_reference_seqs` in the `MutaPipe directory (where the python scripts are stored). The filename should be `UP000005640_9606.fasta` for the script to work properly.
The latest version of this fasta file can be accessed via the [Uniprot ftp client](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/).

**NOTE**: When installing MutaPipe (see [Local Deployment](#local-deployment)) a fasta file of the reference proteome for Homo sapiens will be included (version from December 2021). Please feel free to update to a newer version if one is available (we will also regularly update the reference proteome file on here, but this is currently not done automatically).

