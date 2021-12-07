# MutaPipe - a PDB Pipeline to identify high quality mutant and WT structures for genes of interest
```diff
+TO BE NOTED: We are always working to improve MutaPipe, so any bug reports or suggestions are highly welcome.
```
## Table of Contents
1. [Introduction](#introduction)
2. [Citation](#citation)
3. [Documentation](#documentation)
	* [Minimum Requirements](#minimum-requirements)
	* [Local Deployment](#local-deployment)
	* [Usage](#usage)
	* [Usage Example](#usage-example)
	* [Output](#output)
	* [How to download the reference proteome](#how-to-download-the-reference-proteome)

## Introduction

MutaPipe is a fast and efficient bioinformatics pipeline to screen the Protein Data Bank for genes of interest and retrieve the highest quality protein structure (best resolution) for each unique sequence associated with a gene (for the species Homo sapiens). Additionally, whenever information on variants is available, variants will be annotated using information from ClinVar.

This allows researchers to efficiently screen the PDB for the most suitable template structures for a specific WT or mutant gene/protein which can be used for further in silico analysis like mutagenesis experiments or molecular dynmamics simulations.

![alt text](https://github.com/Utilon/Pipeline_Git/blob/main/debs_files/pipeline_simple_flowchart.jpeg)

Figure 1. Pipeline overview. MutaPipe requires a list of gene names in txt format as input. The pipeline firstly performs a search of the entire PDB database to retrieve information on all protein structures available in the PDB which are associated with the input genes (species: Homo sapiens). In a second step, a BLASTp against the canonical reference sequence for a given input gene (obtained from UNIPROT) is performed for each sequence of every identified protein structure associated with said gene of interest. The BLASTp alignment is used to identify mismatches in the amino acid sequence of the PDB structures relative to the canonical sequence. In a subsequent step, the data is filtered and rearranged to output a table listing the highest quality structure for all available sequences in the PDB which are associated with the input genes. Additionally, if a given amino acid change in an input gene is listed in ClinVar, this information will be added to the output table.


Here I am going to write some more about the pipeline and maybe will also include the table (like in my word file) of the scripts and what they do (needs to be adapted to fit new names etc.)


![alt text](https://github.com/Utilon/Pipeline_Git/blob/main/debs_files/Pipeline_detailled_workflow.jpeg)

Figure 2. Detailled description of the analyses performed by MutaPipe


## Citation

This work has not yet been published (Dec 2021). However, a dataset produced with Mutapipe for all genes available in the ClinVar Database has been presented as a poster at an international research conference in Italy in September 2021 ("THE 1° MASBIC - DISVA ANNUAL SYMPOSIUM - From structure to function: unveiling the role of proteins in health and disease")

## Documentation

### Minimum Requirements

- python >= 3.6 (?)
- RAM: 
- Space required by the installation: 
- Scratch space for usage: 

#### Local Deployment

To obtain MutaPipe please use git to download the most recent development tree:

```bash
git clone https://github.com/Utilon/Pipeline_Git/tree/main/MutaPipe
```

### Usage

In order to run the scripts in the pipeline, a fasta file containing reference proteomes / canonical sequences for the organism in question (=Homo sapiens) will have to be downloaded from UNIPROT and stored in a folder called Uniprot_reference_seqs which should be in the same directory as the scripts. The filename should be UP000005640_9606.fasta for the script to work properly.
The latest version of this fasta file can be accessed via the Uniprot ftp client: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/

When installing MutaPipe (see Local Deployment below) this file will be included (version from December 2021). Do feel free to update with a newer version if one is available (we will also regularly update the reference proteome file on here, but this is currently not done automatically)

IMPORTANT: DNAscan.py is the main script performing the analyses. It must be in the same folder as paths_and_configs.py. Before running DNAscan please modify paths_and_configs.py to match your dependencies deplyment.
IMPORTANT2: All paths in DNAscan end with "/"

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

#### Usage Example

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


### Output

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
           

### How to download the reference proteome

MutaPipe will need a fasta file with all canonical protein sequences of the species of interest.

In order to run the scripts in the pipeline, a fasta file containing reference proteomes / canonical sequences for the organism in question (=Homo sapiens) will have to be downloaded from UNIPROT and stored in a folder called Uniprot_reference_seqs which should be in the same directory as the python scripts. The filename should be UP000005640_9606.fasta for MutaPipe to work properly.
The latest version can be accessed via the Uniprot ftp client: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/

When installing MutaPipe (see Local Deployment below) this file will be included (version from November 2021). Do feel free to update with a newer version if one is available (we will also regularly update the reference proteome file on here, but this is currently not done automatically)

