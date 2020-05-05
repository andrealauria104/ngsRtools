#!/usr/bin/env Rscript
#
# ngsRtools #
#
# Install dependencies ----
#
# 0. Resources ----
CRAN <- "https://cran.rstudio.com/"
#
# 1. CRAN ----
cran_pkgs <- c( 'xlsx', 'ggplot2','ggrepel','reshape2', 'RColorBrewer',
                'circlize', 'plyr', 'ggpubr','ggsci','grid','scales',
                'VennDiagram','Rtsne','rjson','cluster','vegan',
                'dynamicTreeCut', 'gProfileR',"DT")

cran_not_installed <- cran_pkgs[which(!cran_pkgs %in% rownames(installed.packages()))]

if(length(cran_not_installed)>0) {
   for(i in cran_not_installed) {
     install.packages(i, repos=CRAN)	
 }
}
rm(i)
# 2. Bioconductor ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = CRAN)

bioconductor_pkgs <- c('ComplexHeatmap','Biostrings','GenomicRanges','genomation',
                       'BiocGenerics', 'scater', 'scran', 'SC3', 'monocle',
                       'SingleCellExperiment','edgeR','DESeq2', 'sva','clusterProfiler',
                       'org.Mm.eg.db', 'org.Hs.eg.db','methylKit','S4Vectors','IRanges')

bioconductor_not_installed <- bioconductor_pkgs[which(!bioconductor_pkgs %in% rownames(installed.packages()))]

if(length(bioconductor_not_installed)>0) {
   for(i in bioconductor_not_installed) {
     BiocManager::install(i)
 }
}

# 3. Check installed packages ----
cran_not_installed <- cran_pkgs[which(!cran_pkgs %in% rownames(installed.packages()))]
bioconductor_not_installed <- bioconductor_pkgs[which(!bioconductor_pkgs %in% rownames(installed.packages()))]

if(length(cran_not_installed)>0 || length(bioconductor_not_installed)>0) {
  stop(message("[!] Some packages were not successfully installed."))
}
