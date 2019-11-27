# ngsRtools
A Suite of R Packages for NGS-based Epigenomic Data Analysis.

### Description
A Collection of Statistical and Visualization Tools for Next-Generation Sequencing Data Analysis, with a focus on Epigenomics research.<br/> 
<br/>
Includes packages for the analysis of:<br/>

- RNA-seq (bulk transcriptomics) - `RNAseqRtools` <br/>
- scRNA-seq (single-cell transcriptomics) - `scRNAseqRtools` <br/>
- ChIP-seq (genome-wide Protein-DNA binding) - `ChIPseqRtools` <br/>
- BS-seq (genome-wide DNA methylation analysis) - `BSseqRtools` <br/>
- Sequence Data processing and analysis - `sequenceAnalysisRtools` <br/>

### Installation
Clone the repository on your machine:
```
git clone https://github.com/andrealauria104/ngsRtools.git
```
To install the complete suite of packages:
```
cd ngsRtools
./install.sh
```
This will install packages and dependencies. It will also test executability of programs in the `scripts/` folder.<br/>
Alternatively, yuo can install individual packages, for example:
```
R CMD INSTALL RNAseqRtools
```
or from R:
```
devtools::install("RNAseqRtools")
```
To uninstall the complete suite:
```
cd ngsRtools
./uninstall.sh
```
