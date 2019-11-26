# Load and prepare data ----
read_salmon_quant <- function(DATADIR
                              , type = 'gene'
                              , measure = 'counts'
                              , verbose = F)
{
  
  if(type=='gene') {
    
    fl <- list.files(DATADIR
                     , pattern = 'quant.*.sf$'
                     , recursive = T
                     , full.names = T)
    
    nm <- list.files(DATADIR)
    fl <- sapply(nm, function(x)
    {
      tmp <- grep(x, fl, value = T)
      cond <- grepl("genes", tmp)
      if(any(cond)) {
        y <- tmp[which(cond)]
      } else {
        y <- tmp
      }
      return(y)
    }, simplify = T)
  } else if(type=='transcript') {
    
    fl <- list.files(DATADIR
                     , pattern = 'quant.sf$'
                     , recursive = T
                     , full.names = T)
    
  }
  
  get_info <- as.character(system(paste0('cut -f1-3 ', fl[1]), intern = T))
  info <- vector(mode = 'list', length = 3)
  
  for(i in seq_along(info)) {
    info[[i]] <- sapply(strsplit(get_info, "\\t"), "[[",i)[-1]
    names(info)[i] <- sapply(strsplit(get_info, "\\t"), "[[",i)[1]
  }
  info <- do.call(cbind.data.frame, info)
  rm(get_info)
  
  if(measure == 'counts') {
    message("[+] Get counts ...")
    if(verbose) print(fl)
    counts <- lapply(fl, function(i) as.numeric(system(paste0('cut -f5 ', i), intern = T)[-1]))
  } else if(measure == 'TPM') {
    message("[+] Get TPM ...")
    counts <- lapply(fl, function(i) as.numeric(system(paste0('cut -f4 ', i), intern = T)[-1]))
  } else {
    stop(message("[!] Incorret quantification data provided (counts, tpm)."))
  }
  
  sub_idx <- grep('quant.sf', fl)
  if(length(sub_idx)>0) {
    counts[sub_idx] <- lapply(counts[sub_idx], function(i) rep(0,length(counts[[1]])))
  }
 
  names(counts) <- nm
  counts <- do.call(cbind, counts)
  rownames(counts) <- info$Name
  if(any(grepl("NA", colnames(counts)))) {
    counts <- counts[,-which(grepl("NA", colnames(counts)))]
  }
  return(counts)
}

build_expression_matrix <- function(path_quant, path_qc_cells
                                    , pipeline = "hisat"
                                    , measure = "counts"
                                    , rm.ne = T
                                    , ...)
{
  message("[+] Loading expression data, pipeline = ", pipeline)
  get_quant <- function(measure, path_quant, rm.ne, ...)
  {
    if(pipeline == "hisat") {
      m <- loadFunc(type = toupper(measure), countsDIR = path_quant, ...)
    } else if(pipeline == "salmon") {
      m <- read_salmon_quant(DATADIR = path_quant, measure = measure)
    }
    
    colnames(m) <- gsub("_trimmed","",colnames(m))
    m <- m[, qc_cells$sample]
    c.idx <- grep("NA|EMPTY|Minibulk|x", colnames(m), invert = T)
    m <- m[, c.idx]
    if(rm.ne) m <- m[rowSums(m > 0) > 0,]
  }
  
  qc_cells <- readRDS(path_qc_cells)
  rownames(qc_cells) <- gsub("_trimmed","",rownames(qc_cells))
  
  m <- lapply(measure, get_quant, path_quant = path_quant, rm.ne = rm.ne, ...)
  names(m) <- measure
  
  return(list("quant" = m, "info" = qc_cells))
}

# Prepare monocle object ----
prepare_cds_from_sce <- function(sce, pre.normalized = T)
{
  pd <- new("AnnotatedDataFrame", data = as.data.frame(colData(sce)))
  fd <- new("AnnotatedDataFrame", data = as.data.frame(rowData(sce)))
  
  if(pre.normalized) {
    # Use previous normalization
    cds <- monocle::newCellDataSet(normcounts(sce)
                          , phenoData   = pd
                          , featureData = fd
                          , expressionFamily = VGAM::negbinomial.size())
    sizeFactors(cds) <- sizeFactors(sce)
  } else {
    # Use monocle normalization
    cds <- monocle::newCellDataSet(counts(sce)
                          , phenoData   = pd
                          , featureData = fd
                          , expressionFamily = VGAM::negbinomial.size())
    cds <- estimateSizeFactors(cds)
  }
  
  cds <- estimateDispersions(cds)
  
  return(cds)
}
