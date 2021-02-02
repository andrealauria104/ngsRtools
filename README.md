# ngsRtools
A Suite of R Packages for NGS-based Epigenomic Data Analysis.

### Description
A Collection of Statistical and Visualization Tools for Next-Generation Sequencing Data Analysis, with a focus on Epigenomics research.<br/> 
<br/>
It includes packages for the analysis of:<br/>

- RNA-seq (bulk transcriptomics) - `RNAseqRtools` <br/>
- scRNA-seq (single-cell transcriptomics) - `scRNAseqRtools` <br/>
- ChIP-seq (genome-wide Protein-DNA binding) - `ChIPseqRtools` <br/>
- BS-seq (genome-wide DNA methylation analysis) - `BSseqRtools` <br/>
- Sequence Data processing and analysis - `sequenceAnalysisRtools` <br/>

### Prerequisites
The suite is based on a set of CRAN and R/Bioconductor packages. The complete list of 
required packages is reported here for each component of the suite:<br/>

- `RNAseqRtools`: `edgeR`, `DESeq2`, `sva`, `cluster`, `gProfileR`, `clusterProfiler`, `dynamicTreeCut`, `ComplexHeatmap`.<br/>
- `scRNAseqRtools`: `scater`, `scran`, `monocle`.<br/>
- `ChIPseqRtools`: `GenomicRanges`. <br/>
- `BSseqRtools`: `methylKit`. <br/>
- `sequenceAnalysisRtools`: `Biostrings`. <br/>
- `utilsRtools`: `ggplot2`, `reshape2`, `plyr`. <br/>

The installation script will automatically take care of dependencies (see **Installation** section).

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
Alternatively, yuo can install individual packages, for example typing from your R session:
```
devtools::install_github("https://github.com/andrealauria104/ngsRtools", subdir="packages/RNAseqRtools")
```
Or from the location of you cloned repository: 
```
# from command line
R CMD INSTALL RNAseqRtools

# from R
devtools::install("RNAseqRtools")
```
Warning: all packages depend from the `utilsRtools` package. If you choose to install individual packages, install it as you first.<br/> 
To uninstall the complete suite:
```
cd ngsRtools
./uninstall.sh
```
### Using Conda
Clone the repository on your machine:
```
git clone https://github.com/andrealauria104/ngsRtools.git
```
Create conda environment from .yml file with all dependencies: 
```
cd ngsRtools
conda env create -f data/environment.yml
```
Install the complete suite of packages:
```
cd ngsRtools
./install.sh
```
### Docker
The suite can be used in a Docker container with all packages and scripts. Pull the image:
```
docker pull anlauria/ngsrtools
```
To run interactive bash in current directory
```
docker run --rm -it -u `id -u`:`id -g` -v $(pwd):/tmp/ anlauria/ngsrtools:0.0.1 bash
```
To run interactive R console in current directory
```
docker run --rm -it -u `id -u`:`id -g` -v $(pwd):/tmp/ anlauria/ngsrtools:0.0.1 R
```
To run RStudio server in current directory - localhost:8787, username=rstudio
```
docker run --rm -v $(pwd):/home/rstudio -e PASSWORD="ngsrtools" -p 8787:8787 anlauria/ngsrtools:0.0.1
```
