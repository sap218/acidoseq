# acidoseq

Studying Acidobacteria reads | Beta Version (v0.1) | Python v3.5

Author __Samantha C Pendleton__

Please **note**: running through a terminal doesn't seem to operate, however when I use it with Spyder the scripts work.

The **GC** content of the Acidobacteria genomes are consistent with their placements, e.g. species in the same subdivision (above 60\% for group V fragments and roughly 10\% lower for group III fragments) are similar, displaying the diversity within the phylum [1].

The abundance of the subdivisions correlate with pH depends on the subdivisions: 1, 2, 3, 12, 13 have a negative relationship as pH increases, whilst 4, 6, 7, 10, 11, 16, 17, 18, 22, 25 are sparse in low pH and have a positive relationship as pH increases [2].

This package includes assigning unclassified Acidobacteria reads you have collected into subdivisions based on the pH level of the soils. **Note**: not all set-up yet.

## Introduction
[**Kaiju**](http://kaiju.binf.ku.dk) output provides taxon ID and the corredponding sequence, my package outputs the Acidobacteria species alongside annotation, plots, and information on the unclassified reads.

###### Prerequisite
- FASTA format of all the reads.
- Kaiju output after extracting the two columns: sequence ID and NCBI taxIDs.

## Installation

**GitClone**

Note: no command line interface properly set-up as of yet: installations set up through GitClone and I would suggest inputting the Kaiju output and the FASTA of all reads into the same directory.

**Kaiju**

I used the Kaiju output: columns 2 and 3 which included sequence references and the NCBI taxons.

1. Filter the output with only classified labels	`$ awk '$1 == "C"' kaiju.out > kaijuC.out`
2. Cut the columns					`$ cut -f2,3 kaijuC.out > results.txt`
3. Converted the txt to csv (comma-delimted)		`$ sed 's/\s\+/,/g' results.txt > result_seqid_taxon.csv`

## Usage
Run like followed with **Linux**(find how to [run with other operating systems here](https://en.wikibooks.org/wiki/Python_Programming/Creating_Python_Programs).):

`python acidoseq.py`

**Input**

First the script will require information about Acidobacteria data type (all/unclassified), the location of your Kaiju results plus the FASTA file:

`All species or only unclassified ('all' or 'U')?` `'U'`

`Enter your Kaiju Output (edited) file:` `'result_seqid_taxon.csv'`

`Enter the FASTA file of all reads:` `'all.fa'`

**Output**

All outputs will deposit in the **output** file.

The first output includes a FASTA of the reads which are Acidobacteria, either all or only the unclassified depending on your first input. 

`acido-reads_U_2018-08-13_11-58-06.fa`

The terminal will then provide some statistics of the data: Acidobacteria coverage, read-lengths, min/max of AT/GC, plus the means.

Afterwards, you will be prompted for a plotting style (it provides you a selection of styles):

`Insert the style you want:` `ggplot`

For the subdivisions comparison for the GC ratio you will be asked for what kind of plot:

`Enter plot type ('span' or 'line'):` `'span'`

Finally, you will be asked about pH:

`Insert pH of soil:` `'3.23'`


## Map
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
- [x] Improve descriptions and comments
- [ ] Alter code so the input file can be the original Kaiju output
- [ ] Look into command line interface

## Thank you! :seedling:

Don't hesitate to create an issue or make a suggestion!

###### References
[1] Quaiser, A., Ochsenreiter, T., Lanz, C., Schuster, S. C., Treusch, A. H., Eck, J., & Schleper, C. (2003). Acidobacteria form a coherent but highly diverse group within the bacterial domain: evidence from environmental genomics. Molecular microbiology, 50(2), 563-575.

[2] Eichorst, S. A., Breznak, J. A., & Schmidt, T. M. (2007). Isolation and characterization of soil bacteria that define Terriglobus gen. nov., in the phylum Acidobacteria. Applied and environmental microbiology, 73(8), 2708-2717.
