# ngsRtools
A Suite of R Packages for NGS-based Epigenomics Data Analysis.

### Description
A Collection of Statistical and Visualization Tools for Next-Generation Sequencing Data Analysis, with a focus on Epigenomics research.<br/> 
<br/>
Includes packages for the analysis of:<br/>

- RNA-seq (bulk transcriptomics)<br/>
- scRNA-seq (single-cell transcriptomics)<br/>
- ChIP-seq (genome-wide Protein-DNA binding)<br/>
- BS-seq (genome-wide DNA methylation analysis)<br/>

### Installation
Using devtools:

```
install.packages("devtools")
devtools::install_github("https://github.com/andrealauria104/ngsRtools.git")
```

Alternatively, cloning the repo on your machine:

```
git clone https://github.com/andrealauria104/ngsRtools.git
R CMD INSTALL ngsRtools
```
From R:
```
devtools::instal("ngsRtools")
```
