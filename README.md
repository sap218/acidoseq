# acidoseq

Studying Acidobacteria reads | Beta Version (v0.1)

Author __Samantha C Pendleton__

The **GC** content of the Acidobacteria genomes are consistent with their placements, e.g. species in the same subdivision (above 60\% for group V fragments and roughly 10\% lower for group III fragments) are similar, displaying the diversity within the phylum [1].

The abundance of the subdivisions correlate with pH depends on the subdivisions: 1, 2, 3, 12, 13 have a negative relationship as pH increases, whilst 4, 6, 7, 10, 11, 16, 17, 18, 22, 25 are sparse in low pH and have a positive relationship as pH increases [2].

This package includes assigning unclassified Acidobacteria reads you have collected into subdivisions based on the pH level of the soils. **Note**: not all set-up yet.

## Introduction

[**Kaiju**](http://kaiju.binf.ku.dk) output provides taxon ID and the corredponding sequence, my package outputs the Acidobacteria species alongside annotation, plots, and information on the unclassified reads.

###### Prerequisite
- Input Kaiju Output after extracting the two columns: sequence ID and NCBI taxonomy list (`result_seqid_taxon.csv`) - you will need to gain this yourself!
- Input list of NCBI taxons of Acidobacteria (`acido_taxid.csv`), in folder **input**
- Input the FASTA format: you will need a file of all the reads!
- Output FASTA file of all the matched Acidobacteria reads (`acido_reads_2018-07-28_22-28-17.fa`)
- Input this FASTA file to the other scripts

## Installation

###### Pip
**Note**: no installations set up yet, for now you can copy/download the scripts yourself, just make sure you edit the file paths when you run it, or add the correct files into the same directory for ease of use (rather than include the whole file directory, e.g. `\home\samantha\...`) - run like followed with **Linux**:

`python kaiju_taxon_search.py`

Find how to [run with other operating systems here](https://en.wikibooks.org/wiki/Python_Programming/Creating_Python_Programs).

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

This file is what you should input into the other **Python** scripts, e.g. `acgt-comparison.py` and `gc-ratio.py`.

**Example**

`python location-map.py`

You will then be prompted to locate the csv I provided (`latlon.csv`) and the location of your soil sample, which you should input:

`Insert location of csv filename: "latlon.csv"`

`Insert city (case-sensitive): "Aberystwyth"`

You will then gain an image (e.g. `location_soil-ph_2018-08-12_18-16-57.png`) of in that same file directory. **Note**: see how I input quotation marks around the inputs. You can observe an example of a map plot, Liverpool, in the **example** folder.

## Acknowledgements
* **Amanda Clare**, supervisor at Aberystwyth University, [Twitter](https://twitter.com/afcaber) | [GitHub](https://github.com/amandaclare)
* **Sam Nicholls**, postdoc at Birmingham University, [Twitter](https://twitter.com/samstudio8) | [GitHub](https://github.com/SamStudio8)

###### Todo List
- [x] Make available
- [ ] Improve descriptions and comments
- [ ] Alter code so the input file can be the original Kaiju output
- [ ] Look into command line interface

## Thank you! :seedling:

Don't hesitate to create an issue or make a suggestion!

###### References
[1] Quaiser, A., Ochsenreiter, T., Lanz, C., Schuster, S. C., Treusch, A. H., Eck, J., & Schleper, C. (2003). Acidobacteria form a coherent but highly diverse group within the bacterial domain: evidence from environmental genomics. Molecular microbiology, 50(2), 563-575.

[2] Eichorst, S. A., Breznak, J. A., & Schmidt, T. M. (2007). Isolation and characterization of soil bacteria that define Terriglobus gen. nov., in the phylum Acidobacteria. Applied and environmental microbiology, 73(8), 2708-2717.
