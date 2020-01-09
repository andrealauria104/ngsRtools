#!/usr/bin/env Rscript
#
# ngsRtools #
#
# Install dependencies ----
#
# 1. CRAN ----
cran_pkgs <- c( 'xlsx', 'ggplot2','ggrepel','reshape2', 'RColorBrewer',
                'circlize', 'plyr', 'ggpubr','ggsci','grid','scales',
                'VennDiagram','Rtsne','rjson','cluster','vegan',
                'dynamicTreeCut', 'gProfileR')

cran_not_installed <- cran_pkgs[which(!cran_pkgs %in% rownames(installed.packages()))]

if(length(cran_not_installed)>0) lapply(cran_not_installed, install.packages)

# 2. Bioconductor ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioconductor_pkgs <- c('ComplexHeatmap','Biostrings','GenomicRanges','genomation',
                       'BiocGenerics', 'scater', 'scran', 'SC3', 'monocle',
                       'SingleCellExperiment','edgeR','DESeq2', 'sva','clusterProfiler',
                       'org.Mm.eg.db', 'org.Hs.eg.db','methylKit','S4Vectors','IRanges')

bioconductor_not_installed <- bioconductor_pkgs[which(!bioconductor_pkgs %in% rownames(installed.packages()))]

if(length(bioconductor_not_installed)>0) lapply(bioconductor_not_installed, BiocManager::install)