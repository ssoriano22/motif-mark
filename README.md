motif-mark
# Bi625 Motif-Mark

## Program Description

This program visualizes motif sites in DNA or RNA sequences. A single FASTA file - ending in either ".fasta" or ".fa" - containing up to 10 gene sequences of 1000 bases or less is entered as input to the program. Each sequence in this file must have introns identified with lowercase bases and exons idenified with uppercase bases. Also entered as input is a ".txt" file of up to 5 motif sequences consisting of 10 bases or less, where each motif is one line of the file. 

The command line prompt format is as follows:

    motif-mark-oop.py -f Input_Sequences.fasta -m Motifs_of_Interest.txt

The resulting output is a ".png" image with the same prefix as the input FASTA file (i.e. "Figure_1.fasta" -> "Figure_1.png"). This image will contain: 
* A title based on the input FASTA file
* A legend describing how each element (intron, exon, motifs) is depicted in the image
* A separate figure for each gene contained in the FASTA file, with motifs identified by color. All features are drawn to scale for each gene.

## Assignment Requirements (from Bi625 Canvas):

### Minimum Requirements
----------------------
* Well commented, Python3 compatible, object-oriented code, with CLEAR readme.md file
* Public GitHub repo named motif-mark [WARNING! Assignment will not be graded if repo incorrectly named!]
* Script named motif-mark-oop.py [WARNING! Assignment will not be graded if script incorrectly named!]
* Use argparse with options: [WARNING! Assignment will not be graded if argparse options not as requested!]
    * -f: fasta file
    * -m: motifs file
* Output file has same prefix as input file (e.g. Figure_1.fa -> Figure_1.png)
* Input FASTA file (seqs ≤1000 bases) and motifs file (≤10 bases each, one motif per line in a text file)
* Capable of handling:
    * Motifs with ambiguous nucleotides (see https://en.wikipedia.org/wiki/Nucleic_acid_notationLinks to an external site.)
    * Multiple sequences (max 10 in the data you will be provided)
    * Multiple motifs (max 5 in the data you will be provided)
* Consider:
    * How you will handle overlapping motifs
    * How you will denote introns/exons
* All features (motifs, introns, exons) should be to scale
* Output single, well-labeled figure, per FASTA file
* png output
* Key/labeling
* Should be able to be run in the following environment:
    * conda create -n my_pycairo pycairo
    * conda activate my_pycairo
    * [WARNING! Assignment will not be graded if other packages are required!]

### Stretch goals
----------------------
* Staggered drawing of motifs (along the y-axis) to better show position of overlapping motifs
* Transform information in FASTA header to more readable figure title
