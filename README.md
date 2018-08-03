# acidoseq

Studying Acidobacteria reads | Beta Version (v0.1)

Author __Samantha C Pendleton__

[**Kaiju**](http://kaiju.binf.ku.dk) output provides taxon ID and the corredponding sequence, my package outputs the Acidobacteria species alongside annotation, plots, and information on the unclassified reads.

###### Prerequisite
- Input Kaiju Output after extracting the two columns: sequence ID and NCBI taxonomy list (`result_seqid_taxon.csv`) - you will need to gain this yourself!
- Input list of NCBI taxons of Acidobacteria (`acido_taxid.csv`), in folder **input**
- Input the FASTA format: you will need a file of all the reads!
- Output FASTA file of all the matched Acidobacteria reads (`acido_reads_2018-07-28_22-28-17.fa`)
- Input this FASTA file to the other scripts

## Installation

###### Pip
**Note**: for now you can copy the scripts yourself, just make sure you edit the file paths!

###### Files
I used the Kaiju output: columns 2 and 3 which included sequence references and the NCBI taxons.

1. Filter the output with only classified labels	`$ awk '$1 == "C"' kaiju.out > kaijuC.out`
2. Cut the columns					`$ cut -f2,3 kaijuC.out > results.txt`
3. Converted the txt to csv (comma-delimted)		`$ sed 's/\s\+/,/g' results.txt > result_seqid_taxon.csv`

## Usage
Use this file first: `kaiju_taxon_search.py`

###### Input
Result csv file (`result_seqid_taxon.csv`) of the sequence IDs and NCBI taxons from the Kaiju output. 
List of all the NCBI taxonomy numbers and their corresponding Acidobacteria species (`acido_taxid.csv`).
Plus the **FASTA** of all your reads.

**Note**: you will need to include file directory paths!

###### Output
FASTA file: list of the reads which match the Acidobacteria taxons: `acido_reads_2018-07-28_22-28-17.fa`

This file is what you should input into the other **Python** scripts, e.g. `AGCT.py`.

## Acknowledgements
* **Amanda Clare**, supervisor, [Twitter](https://twitter.com/afcaber) | [GitHub](https://github.com/amandaclare)
* **Sam Nicholls**, developer, [Twitter](https://twitter.com/samstudio8) | [GitHub](https://github.com/SamStudio8)

###### Todo List
- [x] Make available
- [ ] Improve descriptions and comments
- [ ] Alter code so the input file can be the original Kaiju output

## Thank you! :seedling:

Don't hesitate to create an issue or make a suggestion!
