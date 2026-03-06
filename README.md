# The Coli Toolkit (CTK): An extension of the modular Yeast Toolkit for use in _E. coli_


This python package contains the code responsible for clustering small DNA fragments in preparation for *de novo* synthesis. The project is also described in the paper: [**The Coli Toolkit (CTK): An extension of the modular Yeast Toolkit to the E. coli chassis**](https://pubs.acs.org/doi/10.1021/acssynbio.5c00489)
 by **Jacob Mejlsted, Erik Kubaczka, Sebastian Wirth, and Heinz Koeppl.


## Install

`pip install DNA-fragment-clustering`

## Python Usage

`from DNA_fragment_clustering import DNA_clustering, DNA_typer`

`DNA_clustering("input.csv", aggressive = False)`

If `aggressive = True` is used, the algorithm will combine single sequences to achieve a higher level of compression, but this may saccrifice synthesizeability due to sequence similarity.



The DNA typer function only uses the path of a .csv file as input:

`DNA_typer("input.csv")`


It is also possible to use the two functions together:

`DNA_clustering(DNA_typer("input.csv"), aggressive = False)`

## Clustering of _de novo_ DNA fragments

The Python function `DNA_clustering` performs clustering and grouping of *de novo* DNA fragments meant for synthesis. From the methods:

>The clustering software uses the Levenshtein similarity matrix to compute the differences between the various fragments that the user wants to synthesize. Using affinity propagation, the software defines clusters with high sequence similarity. From this, groups are made of up to three sequences from distinct clusters to obtain low sequence similarity in the final DNA sequence sent for synthesis. If the aggressive clustering option is selected, groups only containing one sequence are concatenated together to minimize the amount of DNA needed to be synthetized. Following the grouping, the DNA sequences are concatenated and the restriction sites for BsmBI are exchanged to BbsI and BspMI for the second and third occurrences, respectively. The final sequence is then outputted as a .csv file to the same folder as the input file was chosen from.

### Input format

The input.csv files were based on the output format of Benchling.  
The format requires two columns: **Name**, **Sequence**
These are the name of the DNA fragment, and sequence in question, respectively. All other columns will be ignored


## Typing of DNA fragments

The Python function `DNA_typer` adds bases to the 5'- and 3'-ends of the sequence to determine its part type, and to enable entry cloning into pYTK001.

### Input format
The input.csv file requires three columns: **Name**, **Type**, and **Sequence**.
These are the name of the DNA part, , it's type according to the YTK/CTK nomenclature, and the sequence in question, respectively. All other columns will be ignored

## Citation
If you use this code or the data provided here, please cite the corresponding paper. 


## License
The code and the data is available under an MIT License. Please cite the corresponding paper if you use our code and/or data.

## Funding & Acknowledgments
The authors acknowledge Anika Kofod Petersen for her work on the prototype of the _de novo_ synthesis clustering pipeline. 
The work was made possible with the support of a scholarship from the German Academic Exchange Service (DAAD), project number 91877921 to J.M. E.K. was supported by ERC-PoC grant PLATE (101082333). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the funding agencies.
We acknowledge the use of Python and the aforementioned Python packages.
