# CoRINs (Comparator of Residue Interaction Networks)

CoRINs is a web tool for comparing multiple residue interaction networks (RINs) that generates easy to understand results through interactive tables, graphics and downloable text files. CoRINs reports may be useful in the evaluation of protein conformational variation, as an additional tool to validate models from homology modeling, or to assess the putative impacts of a mutation on the protein structure and activity, by identifying its effects on the protein network.

CoRINs does not construct RINs. To construct the Residue Interaction Networks (RINs) from PDB protein structure files, you need to use the [RING 2.0](https://protein.bio.unipd.it/ring/) program or webserver (Piovesan et al. 2016 - [PMID: 27198219](https://pubmed.ncbi.nlm.nih.gov/27198219/)). For now, CoRINs only uses [RING 2.0](https://protein.bio.unipd.it/ring/) output files as input, but we pretend to expand this option soon.

The tool's central core uses Python language with high-performance libraries such as Pandas (https://pandas.pydata.org/), Numpy (https://numpy.org), and Matplotlib (https://matplotlib.org/). For the Web modules, CoRINs uses the frameworks Django (www.djangoproject.com/), Bootstrap (https://getbootstrap.com/), and D3.js (https://d3js.org/). There is also a graph parameter generation module (betweenness and clustering coefficient) developed in the R language (https://www.r-project.org/) and using R package igraph (https://igraph.org), which calculates betweenness weighted and clustering coefficient.

## Authors: 
 
 - Felipe Vieira da Fonseca (Author, Data Scientist)
 - Thiago Dantas Soares (Web Developer, Technical Support)
 - Romildo O. Souza Júnior (Data Scientist, Technical Support)
 - Diego Morais (R Developer)
 - João Paulo Matos Santos Lima (Project supervisor).

Bioinformatics Graduate Program, Instituto Metrópole Digital ([IMD](https://imd.ufrn.br/portal/)), Bioinformatics Multidisciplinary Environment ([BioME](bioinfo.imd.ufrn.br), Universidade Federal do Rio Grande do Norte ([UFRN](ufrn.br)).

### Citation: 

Fonseca, FV, Souza Júnior, RO, Almeida, MVA, Soares, TD, Morais, DAA, Dalmolin, RJS, Lima, JPMS. CoRINs: A tool to compare residue interaction networks from homologous proteins and conformers. BioRxiv. 2020. DOI:[10.1101/2020.06.29.178541](https://www.biorxiv.org/content/10.1101/2020.06.29.178541v1).

[https://www.biorxiv.org/content/10.1101/2020.06.29.178541v1](https://www.biorxiv.org/content/10.1101/2020.06.29.178541v1)

## How does it work?

CoRINs uses as input [RING 2.0](https://protein.bio.unipd.it/ring/) output files, ```nodes.txt```, and ```edges.txt```, or their compressed form (.zip file). The nodes and edges are converted into data frames, separated by protein chains, and pairwisely compared according to the amino acid position, and the different protein chains. The program compares each protein chain, including the one from the same protein.

From the initial screen of the CoRINs, you can submit all the RINs using the Load button, which will open a file manager according to the operating system used. After the load, CoRINs reports the results in interactive graphs and tables according to the type of analysis. 

## Interpreting the results

CoRINs can show some useful information for network analysis in a global way, i.e., a comparison between all loaded and local RINs, comparing RINs chain by chain.

### Standard deviation of degree values from all comparisons:

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico1-standard-deviation-of-degree-of-all-comparisons.png)

This graph displays in the x-axis, the amino acid position, and in the y-axis, the standard deviation of degree values (Dsd) with its result reflected negatively for better visualization of the variation of the degree in that position. Two guiding lines represent the standard deviation values of 0.5 and 1.5 with the following color code (Table).

| Dsd Intervals     | Cor         | 
| :----------------:|:-----------:|
| >= 0.0 and <= 0.5 | Blue        |
| > 0.5 and < 1.5   | Brown       | 
| >= 1.5            | Red         |   
  
Standard deviations around the mean degree value with variation above 0.5 already indicate residues that vary considerably in their number of connections. Standard deviation values greater than 1.5 indicate probable positions of considerable conformational variation.

### Overall differences in interactions (edges):

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico2-Total-differences-in-interactions.png)

The plot of overall differences in interactions consists of the total number of chemical interactions that are not present in the comparison of two RINs. CoRINs only shows the interactions described by the [RING 2.0](https://protein.bio.unipd.it/ring/) software (VDW, HBOND, PIPISTACK, PICATION, and IONIC). The x-axis shows the type of interaction and the y-axis the total number of interactions lost. This graph is generated from the D3.js library and has interactivity with the mouse cursor, changing the color of the box from blue to red.

## Degree values

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico3-High-number-degree.png)
  
For better visualization of degree values, CoRINs generates two plots in "heatmap-like" format. It is one plot for each RIN of the comparison, to bring a full visualization of the degrees of each amino acid in the RIN and guide the user on which residues of higher degrees or with higher differences in their values.

When you put the mouse cursor over the heatmap, the amino acid Id, chain, position, and acronym of the amino acid name, as well as its degree value, is shown. The color scale depends on the value associated with the degree of that position (Image), from the lightest, near white, the absence of amino acid in the position, to a darker color, to red, for the higher degrees.

  <img src="https://raw.githubusercontent.com/d3/d3-scale-chromatic/master/img/OrRd.png" width="100%" height="40" alt="OrRd">
  
On the x-axis you can see the number of the first amino acid position in the line and on the y-axis display a sequence of numbers from 0 to the square root of the total number of amino acids. As an example: to locate the amino acid on position 124, search the line with the closest number, 121. Then add the number below the line, in the case the value 3.

### Highest degree differences

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico4-high-numbers-of-degree-differences.png)

The user can get a quick overlook on which amino acids varied the most in the degree values. To do this, CoRINs generates a plot using the D3.js library with the ten most significant differences in network degrees. The x-axis is the amino acid Id, and the y-axis being the difference value. In the RIN from PDB 1B55 chain A, the leucine at position 168 has a degree value of 2. The same amino acid has a degree value of 6, in the RIN calculated from the PDB 2Z0P. Therefore, the plot displays this difference. 

### Interactive Tables (Differences of Degrees, Differences in Aminoacid, Differences of Nodes, Differences of Edges, Graph parameters)

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico5-differences-of-degrees.png)

To better represent the results, CoRINs shows several interactive tables, implemented with the datatables library (https://datatables.net/). Thus, it is possible to verify which amino acids vary in degree, amino acid exchanges, absence/presence of residues, and the clustering coefficient and the betweenness weighted parameters of each amino acid in the network. All the tables are searchable and customizable.

### Overall RINs Differences

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico6-visualization-of-differences-RIN.png)

Another heatmap-like plot displays an overview of the RIN comparison.
In this graph it is possible to observe for different colors, the positions that are with missing amino acids in the chain, the unique presence of an amino acid in only one of the RINs, exchange of amino acids in a specific position and the number of differences in chemical interactions (Table). If you intend to use this plot, we advise you to mark the option to remove them in the RING 2.0 software, since ligands and heteroatoms tend to have a high number of connections. 

| Description       | Color          | 
| ---------------- |:------------:|
| Absent amino acid in the two RINs | Gray         |
| Number of differences of the residues chemical interactions, at that position    | Scale orange-rose       | 
| Amino acid present only in the RIN 1  | Dark Blue     |   
| Amino acid present only in the RIN 2  | Light Blue    |  
| Position in which amino acid exchange occurred  | Red     | 

## Installation

### Windows systems

Get Docker application at [https://www.docker.com/](https://www.docker.com/)

Open CMD, enter the project folder and execute docker build command:

```
docker-compose build
```

Type the command to run the application:

```
docker-compose up
```

### Linux systems

Install Docker:

```
sudo apt-get install -y docker.io
```

Install Docker-Compose:

```
sudo curl -L "https://github.com/docker/compose/releases/download/1.26.0/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
sudo usermod -aG docker "$USER"
```

Open CMD, enter the project folder and execute docker build command:

```
docker-compose build
```

Type the command to run the application:

```
docker-compose up
```
