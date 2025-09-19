## I.	DESCRIPTION
Integrase on Demand (IOD) is a tool to identify candidate integrases and their attachment (attB) sites in any genome of interest (GOI). IOD uses the integrases and attB sites from genomic islands predicted by TIGER and Islander. The objective of this software is to generate a list of integrases that are putatively able to act on an attachment site found in a GOI for validation of the integrase and use as a genetic engineering tool. This is accomplished by searching for attB sequences in the genome of interest. Integrase on demand returns the attachment sequence, coordinates, and integrase(s) predicted to act on the attB sequence (by TIGER/Islander)

## II.	INSTALLATION

The following files are needed to run integrase on demand:

Two python scripts in bin/ folder:
- runiod.py
- iod.py

Four reference files in lib/ folder:
- reflist.txt
- isles.pkl
- reps.msh
- ints.gff

Note for reviewers: The files needed for integrase on demand to run is available for download at Zenodo.
1. Navigate to the following link https://zenodo.org/records/17088441?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjJlMTVmNzc0LTlmNDItNDFjYS04MjhiLWE0NDhjYTlhYTlmMiIsImRhdGEiOnt9LCJyYW5kb20iOiIwNTJkMzA3YWYzYzYxZWUwZmJkZGY3M2MyMThjNmJlOSJ9.jb-E8lseIiv83XZAF9PaFirM45Z_sw3FIi_PRi_CVPY4obwIoYPiMEu2U2goqO9pKRTKXVPunM2ufUv42mvNhw
2. Scroll down to files and select "Download all"
3. Move downloaded files into the folder: Integrase-On-Demand/lib/

Script to automate dataset download will be available soon!

Dependencies
- python 3.12
- mash 2.3 (https://github.com/marbl/Mash)
- BLAST+ 2.16 (https://anaconda.org/bioconda/blast)

Make a new conda environment from scratch:
```
conda create -n new_conda_environment python=3.12
```
```
conda activate new_conda_environment
```
```
conda install -c bioconda mash
```
```
conda install -c bioconda blast
```
________________________________________________________________________________ 
## III.	USAGE GUIDE

Integrase on demand has a single required argument. Arguments and options can be written in any order. The simplest functional command:
```
python <PATH_TO_IOD>/bin/runiod.py [-options] -in genome.fa
```
### Required Argument
-in \<PATH\> 	File path to genome of interest (GOI)

#### Optional Arguments & Flags
-MODE (Optional, Default MODE="tax") Determines the workflow, MODE = "tax" or "search".

-out \<PATH to output directory\> 	(Optional, Default = “./iod_output/”) Path to file directory to write outputs to. If not specified, a directory is made in the current working directory. If path exists it will be written over.

-yl \<Integer\> 	(Optional, Default= 10) Integer value for length of flanks in attB sites of tyrosine integrases

-sl \<Integer\>	(Optional, Default= 16) Integer value for length of flanks in attB sites of serine integrases or s-core integrases

-threads \<Integer\>	(Optional, Default= 1) Maximum number of searches to run in parallel

-n \<Integer\>	(Optional, Default= 500) Maximum number of attBs to search

-seq \<comma separated list OR FILE\> (Optional) AttB sequences to search when using "search" mode

#### Modes

##### Taxonomic
Taxonomic mode runs the MASH program to determine a rough species estimate, by comparing a simplified "sketch" of the input genome against a database of MASH sketches. This reduces the number of attB sequences to query, thereby decreasing runtime and overall memory usage while still yielding the majority of hits when querying all attBs.

Example command:
```
python <PATH>/bin/runiod.py -in genome.fa -out output_folder_name -tax
```
##### Search
Search mode is bifunctional. The user can use search mode without specifying any attB sequences, as in:
```
python <PATH>/bin/runiod.py -in genome.fa -out output_folder_name -search
```
In this case, the entire database of sequences (182,734) are searched. Alternatively the user can specify a list of attB sequences to search. The list may be provided as a comma separated list in standard input:
```
python <PATH>/bin/runiod.py -in genome.fa -out output_folder_name -search -seq attB1,attB2,...,attBn
```
The list of sequences can also be provided in a file:
```
python <PATH>/bin/runiod.py -in genome.fa -out output_folder_name -search -seq path/to/attB_sequences
```
If a path to a file is given, it must be in one of two formats:

- FASTA format

\>attb1

AGCGGTCGCGtccatacGTGTGCGTG

\>attb2

TTAGCGGTGGatccataGTAAGTAAG

- Each sequence on a new line
  

AGCGGTCGCGtccatacGTGTGCGTG

TTAGCGGTGGatccataGTAAGTAAG

### Use-case examples by mode: 
I.	User wants to develop a system for inserting a genetic element in the genome of a particular organism:
  - Use taxonomic mode to find reference islands in near neighbor genomes. Access final_candidates for list of best attB sites and predicted integrase proteins
    
II. User has a particular integrase protein and wants to search for a known corresponding attachment site in the genome:
- Run integrase on demand in search mode. The "-seq" flag should be followed by the the attB sequence or a file containing attB sequence(s)

III. User has a genome which has not been widely studied and thus gets very few or no attB hits in taxonomic mode:
  - Use search mode without "-seq" flag 
    ________________________________________________________________________________
## IV. OUTPUT FORMAT

 Along with a log file (iod_output/iod.log), the program produces:
    1. iod_output/final_candidates.tsv
    2. iod_output/isles.json
    3. iod_output/iod_dupes.tsv
    4. iod_output/ints.gff

In search mode, the following is also produced:
    5. iod_output/iod_islands.tsv

#### 1. final_candidates.tsv
Tab seperated file with the following columns:

- Contig hit in input genome
- Direction (plus/minus)
- Start-stop coordinates (1-based index)
- attB sequence in input genome
- Example reference island (taken from the highest support island)
    - Integrase
    - IslandID
    - Contig and coordinates of the example reference island
    - Program (TIGER/Islander) source of predicted island
    - Support value
    - Island type (phage, ICE, etc.)

#### 2. isles.json
Reference genomic island information in json format with the following hierarchy:

    attB sequence:
        GTGB-genus__GTDB-species: 
                Integrase,
                IslandID,
                Scaffold and coordinates of genomic island,
                Strand(+/-);attL coordinates;att    coordinates,
                Predicted crossover site coordinates,
                attP,
                Source,
                Support,
                Island type
    

Note, the IslandID contains the first nine digits are GCA accession and the type of gene the island is integrated in the reference. All t(m)RNA genes are represented by a single letter.

#### 3. attb_dupes.tsv
Tab seperated file with the following columns:
- Duplicate category (dupe, ID_block, tandem, or low_supp): representative hit-coordinates
- Contig hit in input genome
- Direction (plus/minus)
- Start-stop coordinates (1-based index)
- attB sequence in input genome
- Each reference island on a tab with the following in a comma separated list:
    - Integrase
    - IslandID
    - Contig and coordinates of the example reference island
    - Program (TIGER/Islander) source of predicted island
    - Support value
    - Island type (phage, ICE, etc.)
 
#### 4. ints.gff
GFF formatted entry fir each integrase associated with each attB sequence in final_candidates.tsv and attb_dupes.tsv. Standard GFF fields with the attributes column containing the following key-value pairs:
  - ID= Unique protein ID (9-digit GenBank accession + "_" + gene number annotated by Prodigal)
  - gca=9-digit GenBank accession
  - name=Shorthand integrase nickname (not unique) based on Pfam domain hitting the protein sequence and the gene number annotated by Prodigal
  - seq=Amino acid sequence

#### 5. occupied.tsv
Tab separated file with the following columns:
- Contig and coordinates
- attL and attR start-stop (1-based index)
- Number of reference islands which mapped to this locus
- Each reference island on a tab with the following in a comma separated list:
    - Integrase
    - IslandID
    - Contig and coordinates of the example reference island
    - Program (TIGER/Islander) source of predicted island
    - Support value
    - Island type (phage, ICE, etc.)

## V.	COPYRIGHTS
Copyright 2025 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
