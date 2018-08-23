# acidoseq

Studying Acidobacteria reads v1 | Python v3.5 

Author __Samantha C Pendleton__, Data Science MSc student at Aberystwyth University, [Twitter](https://twitter.com/sap218) | [GitHub](https://github.com/sap218)

The **GC** content of the Acidobacteria genomes are consistent with their placements, e.g. species in the same subdivision (above 60\% for group V fragments and roughly 10\% lower for group III fragments) are similar, displaying the diversity within the phylum [1].

The abundance of the subdivisions correlate with pH depends on the subdivisions: 1, 2, 3, 12, 13 have a negative relationship as pH increases, whilst 4, 6, 7, 10, 11, 16, 17, 18, 22, 25 are sparse in low pH and have a positive relationship as pH increases [2].

This package includes assigning unclassified Acidobacteria reads you have collected into subdivisions based on the pH level of the soils. **Note**: not all set-up yet.

## Introduction
[**Kaiju**](http://kaiju.binf.ku.dk) output provides taxon ID and the corredponding sequence, my package outputs the Acidobacteria species alongside annotation, plots, and information on the unclassified reads.

###### Prerequisite
- FASTA format of all the reads.
- Kaiju output after extracting the two columns: sequence ID and NCBI taxIDs.

I would store these in the working directory, but not specifically in the input folder.

## Installation

**GitClone**

Note: no command line interface properly set-up as of yet: installations set up through GitClone and I would suggest inputting the Kaiju output and the FASTA of all reads into the same directory.

**Kaiju**

I used the Kaiju output: columns 2 and 3 which included sequence references and the NCBI taxons.

1. Filter the output with only classified labels	`$ awk '$1 == "C"' kaiju.out > kaijuC.out`
2. Cut the columns					`$ cut -f2,3 kaijuC.out > results.txt`
3. Converted the txt to csv (comma-delimted)		`$ sed 's/\s\+/,/g' results.txt > result_seqid_taxon.csv`

## Usage
Run like followed with **Linux** (find how to [run with other operating systems here](https://en.wikibooks.org/wiki/Python_Programming/Creating_Python_Programs).):

```
$ python3 acidoseq.py --help
Usage: acidoseq.py [OPTIONS]

Options:
  --taxdumptype TEXT  Study "ALL" or only unclassified "U"?
  --kaijufile TEXT    Place edited Kaiju (csv) in directory for ease.
  --fastapath TEXT    Place FASTA in directory for ease.
  --style TEXT        ['seaborn-bright', 'seaborn-poster', 'seaborn-white',
                      'bmh', 'seaborn-darkgrid', 'seaborn-pastel',
                      'grayscale', '_classic_test', 'ggplot', 'seaborn-
                      whitegrid', 'seaborn-dark', 'seaborn-muted', 'seaborn-
                      colorblind', 'seaborn-ticks', 'Solarize_Light2',
                      'seaborn-notebook', 'dark_background', 'fast',
                      'seaborn', 'fivethirtyeight', 'seaborn-paper', 'seaborn-
                      dark-palette', 'seaborn-talk', 'classic', 'seaborn-
                      deep']
  --plottype TEXT     "span" range of GC means OR "line" average mean GC
  --ph TEXT           pH of soil, use map script for assistance.
  --help              Show this message and exit.
```

**Example**

`$ python3 acidoseq.py --taxdumptype ALL --kaijufile result_seqid_taxon.csv --fastapath all.fa --style ggplot --plottype span --ph 4.92`

## Map
Please **note**: due to the fact that the Earth is spherical and maps are 2-dimensional, there will be some distortion when plotting locations.

Script is set-up so input city can be any case (e.g. `liVERpool`). You will then gain an image in the output folder, which you can observe an example already. 

`python3 location-map.py`

`Insert city (e.g. Aberystwyth):` `Birmingham`

## Acknowledgements
* **Amanda Clare**, supervisor at Aberystwyth University, [Twitter](https://twitter.com/afcaber) | [GitHub](https://github.com/amandaclare)
* **Sam Nicholls**, postdoc at Birmingham University, [Twitter](https://twitter.com/samstudio8) | [GitHub](https://github.com/SamStudio8)

###### Todo List
- [x] Make available
- [x] Improve descriptions and comments
- [x] Look into command line interface
- [ ] Fix code to output unclassified subdivisions based on pH
- [ ] Alter code so the input file can be the original Kaiju output

## Thank you! :seedling:

Don't hesitate to create an issue or make a suggestion!

###### References
[1] Quaiser, A., Ochsenreiter, T., Lanz, C., Schuster, S. C., Treusch, A. H., Eck, J., & Schleper, C. (2003). Acidobacteria form a coherent but highly diverse group within the bacterial domain: evidence from environmental genomics. Molecular microbiology, 50(2), 563-575.

[2] Eichorst, S. A., Breznak, J. A., & Schmidt, T. M. (2007). Isolation and characterization of soil bacteria that define Terriglobus gen. nov., in the phylum Acidobacteria. Applied and environmental microbiology, 73(8), 2708-2717.
