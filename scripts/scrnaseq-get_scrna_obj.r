#!/usr/bin/env Rscript
# # # # # # # # # # # # # # # # # # # # # # # #
#                                             #
#  Single cell RNA-sequencing data analysis   #
#                                             #
# # # # # # # # # # # # # # # # # # # # # # # #
#
# scrnaseq-get_scrna_obj.r
#
# Create data structure object for downstream analysis
# 
# Supported worflows: Bioconductor (OSCA), Monocle3, Seurat
#
# 0. Resources ----
suppressWarnings(suppressMessages(library(docopt)))

'Single cell RNA-sequencing data analysis
 
 scrnaseq-get_scrna_obj.r
 
 Create data structure object for downstream analysis
 Supported worflows: Bioconductor (OSCA), Monocle3, Seurat

Usage:
   scrnaseq-get_scrna_obj.r [-c <counts> -m <metadata> -w <workflow> -o <outdir> -q <qcmatrix> -g <gencode> -a <analysis>]
              
Options:
   -c, --counts Path to raw counts (comma separated if > 1) (required). 
   -m, --metadata Path to experiment metadata (comma separated if > 1) (required). 
   -g, --gencode Path to gene annotation info (gencode) (optional).
   -q, --qcmatrix Path to QC matrix or multiQC output directory (optional).
   -w, --workflow Workflow used for data processing (bioconductor/monocle/seurat) [default: bioconductor].
   -p, --protocol Protocol used for data generation (smartseq2/sciseq/10x) [default: smartseq2].
   -a, --analysis Project title [default: scRNA-seq].
   -o, --outdir Path to output directory [default: .].

Author:
   Andrea Lauria' -> doc

opts <- docopt(doc)
# required arguments ---
required_args <- opts[1:2]
if(any(sapply(required_args, is.null))) {
  missing_idx <- sapply(required_args, is.null)
  missing_args <- gsub("--"," ",names(required_args)[missing_idx])
  message("\n[!] Missing required arguments: ",missing_args,"\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

# check arguments ---
if(!opts$workflow%in%c("bioconductor","monocle","seurat")) {
  message("\n[!] Invalid workflow\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

if(!opts$protocol%in%c("smartseq2","sciseq","10x")) {
  message("\n[!] Invalid protocol\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}
# resources ----
suppressWarnings(suppressMessages(library(plyr)))

# 1. Parse arguments ----
opts_idx <- c("counts","metadata")
if(!is.null(opts$gencode))  opts_idx <- c(opts_idx,"gencode")
if(!is.null(opts$qcmatrix)) opts_idx <- c(opts_idx,"qcmatrix")

opts[opts_idx] <- lapply(opts[opts_idx], function(x) 
  {
  x <- as.character(x)
  if(any(grepl("\\,",x))) {
    x <- unlist(strsplit(x,"\\,"))
  }
  return(x)
  }
)

# 2. Read data ----
if(as.character(opts$protocol)%in%c("smartseq2","sciseq")) {
  # parse protocols/pipelines generating tabular data 
  read_data_list <- list()
  for(i in opts_idx) {
    if(length(opts[[i]])>1) {
      # integrate multiple runs/datasets
      read_data_list_tmp <- lapply(opts[[i]], function(x) {
        message(" -- reading: ",x)
        read.delim(x
                   , stringsAsFactors = F
                   , header = T
                   , sep = ifelse(grepl("\\.csv$",x),",","\t"))
        
      })
      if(i=="counts" || (i=="gencode" && length(unique(opts[[i]]))!=1)) {
        read_data_list[[i]] <- Reduce(
          function(x, y, ...) merge(x, y, all = TRUE, by = 1),
          read_data_list_tmp
        )
      } else if(i=="metadata" || i=="qcmatrix") {
        read_data_list[[i]] <- do.call(rbind.data.frame, read_data_list_tmp)
      } else if(i=="gencode" && length(unique(opts[[i]]))==1){
        # N.B.: gencode annotation should be the same!
        read_data_list[[i]] <- read_data_list_tmp
      }
    } else {
      message(" -- reading: ",opts[[i]])
      read_data_list[[i]] <- read.delim(opts[[i]]
                                        , stringsAsFactors = F
                                        , header = T
                                        , sep = ifelse(grepl("\\.csv$",opts[[i]]),",","\t"))  
    }
  }
  
  rownames(read_data_list$counts) <- read_data_list$counts[,1]
  read_data_list$counts <- read_data_list$counts[,-1]
  if(!is.null(opts$qcmatrix)) read_data_list$metadata <- merge(read_data_list$metadata, read_data_list$qcmatrix, by=1, all = T)
  rownames(read_data_list$metadata) <- read_data_list$metadata[,1]
  read_data_list$metadata <- read_data_list$metadata[colnames(read_data_list$counts),]
  
  if(!is.null(opts$gencode)) {
    read_data_list$gene_info <- read.delim(as.character(opts$gencode), stringsAsFactors = F, header = T)
    read_data_list$gene_info <- ddply(read_data_list$gene_info, .(gene_name),summarize
                                      , gene_id=paste0(gene_id, collapse = ";")
                                      , gene_type=paste0(gene_type, collapse = ";"))
    rownames(read_data_list$gene_info) <- read_data_list$gene_info$gene_name
    read_data_list$gene_info <- read_data_list$gene_info[rownames(read_data_list$counts),]
  }
  # checks ---
  if(!all(read_data_list$metadata$Sample%in%colnames(read_data_list$counts))) {
    message("\n[!] Colnames in count matrix do not match metadata samples!\n")
    message(doc)
    quit(save = "no", status = 0, runLast = TRUE)
  }
}

# 3. Prepare scrnaseq object ----
if(as.character(opts$workflow)=="bioconductor") {
  suppressWarnings(suppressMessages(library(SingleCellExperiment)))
  
  if(as.character(opts$protocol%in%c("smartseq2","sciseq"))) {
    # Smartseq2 / Sci-RNA-seq data
   if(!is.null(read_data_list$gene_info)) {
     sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(read_data_list$counts))
                                                       , colData = read_data_list$metadata
                                                       , rowData = read_data_list$gene_info)
   } else {
     sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(read_data_list$counts))
                                                       , colData = read_data_list$metadata)
   }
  } else {
    # read 10x data
    suppressWarnings(suppressMessages(library(DropletUtils)))
    sce <- read10xCounts(as.character(opts$counts))
  }
  outfile <- paste0(as.character(opts$outdir),"/sce.rds")
  message(" -- saving object to: ", outfile)
  saveRDS(sce, file = outfile)

} else if(as.character(opts$workflow)=="monocle") {
  suppressWarnings(suppressMessages(library(monocle3)))
  
  if(as.character(opts$protocol%in%c("smartseq2","sciseq"))) {
    # Smartseq2 / Sci-RNA-seq data
    if(is.null(read_data_list$gene_info)) {
      read_data_list$gene_info <- data.frame("gene_name"=rownames(read_data_list$counts), row.names = rownames(read_data_list$counts))
    }
    cds <- new_cell_data_set(as.matrix(read_data_list$counts),
                             cell_metadata = read_data_list$metadata,
                             gene_metadata = read_data_list$gene_info)
  } else {
    # read 10x data
    cds <- load_cellranger_data(as.character(opts$counts))
  }
  outfile <- paste0(as.character(opts$outdir),"/cds.rds")
  message(" -- saving object to: ", outfile)
  saveRDS(cds, file = outfile)
  
} else if(as.character(opts$workflow)=="seurat") {
  suppressWarnings(suppressMessages(library(Seurat)))
  
  if(as.character(opts$protocol%in%c("smartseq2","sciseq"))) {
    # Smartseq2 / Sci-RNA-seq data
    seurat_obj <- CreateSeuratObject(counts = as.matrix(read_data_list$counts)
                                     , meta.data = read_data_list$metadata
                                     , project = as.character(opts$analysis)
                                     , min.cells = 0
                                     , min.features = 0)
    
    if(!is.null(read_data_list$gene_info)) seurat_obj[["RNA"]] <- AddMetaData(seurat_obj[["RNA"]], read_data_list$gene_info)
    
  } else {
    # read 10x data
    scrnaseq_data <- Seurat::Read10X(data.dir  = as.character(opts$counts), project = as.character(opts$analysis))
    seurat_obj    <- CreateSeuratObject(counts = scrnaseq_data)
    
  }
  outfile <- paste0(as.character(opts$outdir),"/seurat_obj.rds")
  message(" -- saving object to: ", outfile)
  saveRDS(seurat_obj, file = outfile)
}
