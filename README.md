## I.	DESCRIPTION
Integrase on Demand (IOD) is a tool to identify candidate integrases and their attachment (attB) sites in any genome of interest (GOI). IOD uses the integrases and attB sites from genomic islands predicted by TIGER and Islander. The objective of this software is to generate a list of integrases that are putatively able to act on an attachment site found in a GOI for validation of the integrase and use as a genetic engineering tool. This is accomplished by searching for attB sequences in the genome of interest. Integrase on demand returns the attachment sequence, coordinates, and integrase(s) predicted to act on the attB sequence (by TIGER/Islander)

## II.	INSTALLATION

The following files are needed to run integrase on demand:

- runiod.py
- iod.py
- isles.pkl
- reps.msh
- reflist.txt

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

Integrase on demand has two required arguments. It can be run by:
```
python <PATH_TO_IOD>/bin/runiod.py -in genome.fa -MODE [-options]
```
#### Arguments
-in \<PATH to genome FASTA\> 	File path to genome of interest (GOI) 
-MODE	\<tax or search\>   Determines the workflow

#### Flags
-out \<PATH to output directory\> 	(Optional, Default = “./iod_output/”) Path to file directory to write outputs to. If not specified, a directory is made in the current working directory. If path exists it will be written over. 	
-yl \<Integer\> 	(Optional, Default =10) Integer value for length of flanks in attB sites of tyrosine integrases
-sl \<Integer\>	(Optional, Default =16) Integer value for length of flanks in attB sites of serine integrases or s-core integrases
-cpus \<Integer\>	(Optional, Default =1) Maximum number of searches to run in parallel
-n \<Integer\>	(Optional, Default =500) Maximum number of attBs to search
-seq \<comma separated list OR FILE\> (Optional) AttB sequences to search

#### Modes

##### Taxonomic
Taxonomic mode runs the MASH program to determine a rough species estimate, by comparing a simplified "sketch" of the input genome against a database of MASH sketches. This reduces the number of attB sequences to query, thereby decreasing runtime and overall memory usage while still yielding the majority of hits when querying all attBs.

Example command:
```
python <PATH>/bin/runiod.py -in genome.fa -out output_folder_name -tax
```
##### Search
Search mode returns for all specified attachment site(s) of interest in the genome.

Example for specifying attBs in standard input:
```
python <PATH>/bin/runiod.py -in genome.fa -out output_folder_name -search -seq attB1,attB2,...,attBn
```
Example command for search mode for fasta file containing query attBs:
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
    
II.	User has a particular integrase protein and wants to search for corresponding attachment site in the genome:
- Run integrase on demand in search mode. The "-seq" flag should be followed by the path to the attB sequences (format can be each sequence on a newline or FASTA)
    ________________________________________________________________________________
## IV. OUTPUT FORMAT

 The program produces:
    1. iod_output/iod_out.tsv
    2. iod_output/isles.json
    3. iod_output/iod_dupes.tsv

In search mode, the following is also produced:
    4. iod_output/iod_islands.tsv

#### 1. iod_out.tsv
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

#### 3.iod_dupes.tsv
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

#### 4.iod_islands.tsv
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
