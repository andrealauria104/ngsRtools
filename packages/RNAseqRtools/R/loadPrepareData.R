# Load and Prepare Data ====
loadCounts <- function(countsDIR, pattern=NULL) {
  
  f  <- list.files(countsDIR        , 
                   pattern    = pattern,
                   full.names = T)
  
  counts <- do.call(cbind, lapply(f, read.delim, header=F, row.names=(1)))
  names(counts) <- gsub(".txt", "", basename(f))
  
  return(counts)
  
  
} 

loadTPM <- function(countsDIR) {
  
  tpmData <- list.files(countsDIR, 
                        pattern    = "TPM", 
                        full.names = T)
  
  tpm <- read.csv2(tpmData, 
                   stringsAsFactors = F,
                   header = T)
  
  x <- lapply(tpm[,-1], as.numeric)
  x <- as.matrix(as.data.frame(x))
  rownames(x) <- tpm$Gene
  
  return(x)
}

loadFunc <- function(type="TPM", ...) {
  
  if( type=="COUNTS" ) {
    return(loadCounts(...)) 
  } else if ( type=="TPM" ) {
    return(loadTPM(...))
  }
}

# still here for compatibility reason, 
# will be deprecated as soon as possible
process_rnaseq_edger <- function(m
                                 , group  = NULL
                                 , reference = NULL # Control group
                                 , filter = T
                                 , filter.expr.th = 1
                                 , filter.sample.th = 2
                                 , normalize.using = "cpm"
                                 , norm.fact.method = "TMM"
                                 , normalized.lib.sizes = TRUE
                                 , tlen = NULL
                                 , ...)
{
  
  message(" -- Pre-process RNA-seq data using edgeR")
  
  if(is.null(group)) {
    group <- as.factor(colnames(m))
  } else {
    m <- m[,names(group)]
    if (!is.factor(group))   group <- as.factor(group)
    if (!is.null(reference)) group <- relevel(group, reference)
  }
  
  message(" -- Condition: ", paste0(levels(group), collapse = "-"))
  y <- edgeR::DGEList(counts=m, genes=rownames(m), group = group)
  message(" -- Normalization factors, method = ", norm.fact.method)
  y <- edgeR::calcNormFactors(y, method =  norm.fact.method, ...)
  if(!normalized.lib.sizes) warning("Normalization factors are ignored in RPKM/CPM calculation. Avoid between samples comparisons!")
  # Clean environment
  rm(m)
  gc(verbose = F)
  
  if(filter) {
    message(" -- Filtering lowly expressed genes")
    message("    -- Threshods:")
    message("     * ",normalize.using," = ", filter.expr.th)
    message("     * samples = " , filter.sample.th)
    
    if(filter.sample.th>=1) {
      # Number of samples
      fs <- filter.sample.th
    } else {
      # Percentage of samples
      fs <- floor(ncol(y)*filter.sample.th)
    }
    if(normalize.using == "cpm") {
      keep <- rowSums(edgeR::cpm(y, normalized.lib.sizes = normalized.lib.sizes)>filter.expr.th) >= fs  
    } else if (normalize.using == "rpkm") {
      keep <- rowSums(edgeR::rpkm(y, gene.length = tlen[rownames(y)], normalized.lib.sizes = normalized.lib.sizes)>filter.expr.th) >= fs  
    } else {
      stop(message("[!] Incorrect normalization method (cpm, rpkm)."))
    }
    
    y <- y[keep, , keep.lib.sizes=FALSE]
  }
  
  if(normalize.using == "cpm") {
    y$CPM    <- edgeR::cpm(y, normalized.lib.sizes = normalized.lib.sizes)
    y$logCPM <- edgeR::cpm(y, normalized.lib.sizes = normalized.lib.sizes, log = T)
  } else if (normalize.using == "rpkm") {
    if(is.null(tlen)) {
      stop(message("[!] Please, provide transcript lenghts for RPKM normalization."))
    }
    y$RPKM    <- edgeR::rpkm(y, gene.length = tlen[rownames(y)], normalized.lib.sizes = normalized.lib.sizes)
    y$logRPKM <- edgeR::rpkm(y, gene.length = tlen[rownames(y)], normalized.lib.sizes = normalized.lib.sizes, log = T)
  }
  
  return(y)
}

processRNAseqEdgeR <- function(m, experimental_info = NULL
                               , gene_info = NULL
                               , group  = NULL
                               , reference = NULL # Control group
                               , filter = T
                               , filter.expr.th = 1
                               , filter.sample.th = 2
                               , normalize.using = "cpm"
                               , norm.fact.method = "TMM"
                               , normalized.lib.sizes = TRUE
                               , tlen = NULL
                               , ...)
{
  
  message(" -- Pre-process RNA-seq data using edgeR")
  
  if(is.null(experimental_info)) {
    if(is.null(group)) {
      group <- as.factor(colnames(m))
      if(sum(table(levels(group))>1)==0) {
        stop(message("[!] Invalid experimental groups (< 2 replicates per condition)."))
      }
    } else {
      if(!is.null(names(group)))  m <- m[,names(group)]
      if(!is.factor(group))   group <- as.factor(group)
      if(!is.null(reference)) group <- relevel(group, reference)
    }
  } else {
    # group as single or combination of columns in experimental info
    # provide as integer, character
    if(!is.null(group)) {
     if(length(group)==1) {
       group <- experimental_info[,group]  
     } else {
       group <- apply(experimental_info[,group], 1, paste0, collapse = "_")
     }
      if(!is.factor(group))   group <- as.factor(group)
      if(!is.null(reference)) group <- relevel(group, reference)
    }
  }
  
  if(is.null(gene_info)) gene_info <- rownames(m)
  
  y <- edgeR::DGEList(counts    = m
                      , genes   = gene_info
                      , samples = experimental_info
                      , group   = group)
  
  message(" -- Condition: ", paste0(levels(y$samples$group), collapse = "-"))
  
  message(" -- Normalization factors, method = ", norm.fact.method)
  y <- edgeR::calcNormFactors(y, method =  norm.fact.method, ...)
  if(!normalized.lib.sizes) warning("Normalization factors are ignored in RPKM/CPM calculation. Avoid between samples comparisons!")
  # Clean environment
  rm(m)
  gc(verbose = F)
  
  if(filter) {
    message(" -- Filtering lowly expressed genes")
    message("    -- Threshods:")
    message("     * CPM = "     , filter.expr.th)
    message("     * Samples = " , filter.sample.th)
    
    if(filter.sample.th>=1) {
      # Number of samples
      fs <- filter.sample.th
    } else {
      # Percentage of samples
      fs <- floor(ncol(m)*filter.sample.th)
    }
    if(normalize.using == "cpm") {
      keep <- rowSums(edgeR::cpm(y, normalized.lib.sizes = normalized.lib.sizes)>filter.expr.th) >= fs  
    } else if (normalize.using == "rpkm") {
      if(is.null(tlen)) {
        stop(message("[!] Please, provide transcript lenghts for RPKM normalization."))
      }
      keep <- rowSums(edgeR::rpkm(y, normalized.lib.sizes = normalized.lib.sizes, gene.length = tlen[rownames(y)])>filter.expr.th) >= fs  
    } else {
      stop(message("[!] Incorrect normalization method (cpm, rpkm)."))
    }
    
    y <- y[keep, , keep.lib.sizes=FALSE]
  }
  
  if(normalize.using == "cpm") {
    y$CPM    <- edgeR::cpm(y, normalized.lib.sizes = normalized.lib.sizes)
    y$logCPM <- edgeR::cpm(y, normalized.lib.sizes = normalized.lib.sizes, log = T)
  } else if (normalize.using == "rpkm") {
    if(is.null(tlen)) {
      stop(message("[!] Please, provide transcript lenghts for RPKM normalization."))
    }
    y$RPKM    <- edgeR::rpkm(y, gene.length = tlen[rownames(y)], normalized.lib.sizes = normalized.lib.sizes)
    y$logRPKM <- edgeR::rpkm(y, gene.length = tlen[rownames(y)], normalized.lib.sizes = normalized.lib.sizes, log = T)
  }
  
  return(y)
}