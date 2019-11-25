# Tools for RNA-seq data analysis

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
                                 , tlen = NULL)
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
  y <- edgeR::calcNormFactors(y)
  
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
      keep <- rowSums(edgeR::cpm(y)>filter.expr.th) >= fs  
    } else if (normalize.using == "rpkm") {
      keep <- rowSums(edgeR::rpkm(y, gene.length = tlen[rownames(y)])>filter.expr.th) >= fs  
    } else {
      stop(message("[!] Incorrect normalization method (cpm, rpkm)."))
    }
    
    y <- y[keep, , keep.lib.sizes=FALSE]
  }
  
  if(normalize.using == "cpm") {
    y$CPM    <- edgeR::cpm(y)
    y$logCPM <- edgeR::cpm(y, log = T)
  } else if (normalize.using == "rpkm") {
    if(is.null(tlen)) {
      stop(message("[!] Please, provide transcript lenghts for RPKM normalization."))
    }
    y$RPKM    <- edgeR::rpkm(y, gene.length = tlen[rownames(y)], normalized.lib.sizes = T)
    y$logRPKM <- edgeR::rpkm(y, gene.length = tlen[rownames(y)], normalized.lib.sizes = T, log = T)
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
                               , tlen = NULL)
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
  }
  
  if(is.null(gene_info)) gene_info <- rownames(m)
  
  y <- edgeR::DGEList(counts    = m
                      , genes   = gene_info
                      , samples = experimental_info
                      , group   = group)
  
  message(" -- Condition: ", paste0(levels(y$samples$group), collapse = "-"))
  
  y <- edgeR::calcNormFactors(y)
  
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
      keep <- rowSums(edgeR::cpm(y)>filter.expr.th) >= fs  
    } else if (normalize.using == "rpkm") {
      if(is.null(tlen)) {
        stop(message("[!] Please, provide transcript lenghts for RPKM normalization."))
      }
      keep <- rowSums(edgeR::rpkm(y, gene.length = tlen[rownames(y)])>filter.expr.th) >= fs  
    } else {
      stop(message("[!] Incorrect normalization method (cpm, rpkm)."))
    }
    
    y <- y[keep, , keep.lib.sizes=FALSE]
  }
  
  if(normalize.using == "cpm") {
    y$CPM    <- edgeR::cpm(y)
    y$logCPM <- edgeR::cpm(y, log = T)
  } else if (normalize.using == "rpkm") {
    if(is.null(tlen)) {
      stop(message("[!] Please, provide transcript lenghts for RPKM normalization."))
    }
    y$RPKM    <- edgeR::rpkm(y, gene.length = tlen[rownames(y)], normalized.lib.sizes = T)
    y$logRPKM <- edgeR::rpkm(y, gene.length = tlen[rownames(y)], normalized.lib.sizes = T, log = T)
  }
  
  return(y)
}
# Quality Controls ====

GrepStats <- function(f, pipeline="Hisat2") {
  
  # Pipeline Hisat2
  if( pipeline=="Hisat2" ) {
    grp <- paste0("grep 'reads\\|alignment\\|aligned' ",f,"  | cut -d ' ' -f1,5")  
  } else {
    stop(message("[-] Please, provide a valid RNA-seq pipeline."))
  }
  
}

AlignStat <- function(STATDIR, ...) {
  
  fl <- list.files(STATDIR
                   , pattern = "txt"
                   , full.names = T)
  
  grps <- vapply(fl, GrepStats, FUN.VALUE = character(1), ...)
  
  tmp <- lapply(grps, system, intern=T)
  names(tmp) <- basename(fl)
  tmp <- do.call(rbind, tmp)
  
  nm <- rownames(tmp)
  
  tmp <- gsub("%| ","",tmp)
  tmp <- as.data.frame(apply(tmp, 2, as.numeric))
  colnames(tmp) <- c("library_size","unmapped","uniquely_map","multi_map","mapping_rate")
  tmp$sample <- gsub("\\.txt","",nm)
  
  return(tmp)
  
}

# Normalization ====
countToTpm <- function(counts, effLen)
{
  effLen <- effLen[rownames(counts)]
  # Transcripts Per Million (TPM)
  tpms <- apply(counts, 2,
                function(x) {
                  rate <- log(x) - log(effLen)
                  denom <- log(sum(exp(rate)))
                  c <- exp(rate - denom + log(1e6))
                  return(c)
                }
  )
  return(tpms)
}

countToFpkm <- function(counts, effLen, pcg)
{
  # Fragments Per Kilobases per Million (FPKM) 
  effLen <- effLen[rownames(counts)]
  pcg    <- intersect(pcg, rownames(counts))
  
  fpkms <- apply(counts, 2,
                 function(x) {
                   N <- sum(x[pcg])
                   c <- exp( log(x) + log(1e9) - log(effLen) - log(N) )
                   return(c)
                 }
  )
  return(fpkms)
}

# Principal Component Analysis ====
plotPCA <- function(x, experimental_info = NULL
                    , labels   = T
                    , col_by   = NULL
                    , shape_by = NULL
                    , pal      = NULL
                    , scree_plot = F
                    , scree_plot_type = "explained_variance"
                    , point_size = 2
                    , max.pcs = NULL
                    , dim_1 = "PC1"
                    , dim_2 = "PC2"
                    , default.title = "RNA-seq PCA") {
  
  rna <- t(x)
  rna <- rna[,colSums(rna)>1]
  RNA_pca <- prcomp(rna, scale. = T, center = T)
  pca_summary <- summary(RNA_pca)$importance
  
  if((is.null(max.pcs) || max.pcs > ncol(RNA_pca$x))) {
    max.pcs <- ncol(RNA_pca$x)
  }
  
  pca_plot <- as.data.frame(RNA_pca$x[,paste0("PC",1:max.pcs)])
  
  if( !is.null(experimental_info) ) {
    if(is.character(experimental_info)) {
      # for compatibility with previous version
      experimental_info <- reshape2::melt(experimental_info)
      if(!is.null(col_by)) {
        colnames(experimental_info)[1] <- col_by 
      }
    }
    if(is.null(col_by)) {
      col_by <- colnames(experimental_info)[1]
    }
    pca_plot[[col_by]] <- experimental_info[rownames(pca_plot), col_by]
    if(!is.null(shape_by)) {
      pca_plot[[shape_by]] <- experimental_info[rownames(pca_plot), shape_by]
    }
  } else {
    pca_plot[[col_by]] <- rownames(pca_plot) 
  }
  
  if(is.null(pal)) {
    pal <- ggsci::pal_d3()
    pal <- pal(length(unique(pca_plot[,col_by])))
  }
  pca_plot$repel_col_by <- pca_plot[,col_by]
  
  p <- ggplot(pca_plot, aes_string(x=dim_1, y=dim_2, col=col_by, shape=shape_by)) + geom_point(size=point_size) +
    xlab(paste0(dim_1," (",round(pca_summary[2,grep(paste0(dim_1,"$"), colnames(pca_summary))]*100,1),"%)")) +
    ylab(paste0(dim_2," (",round(pca_summary[2,grep(paste0(dim_2,"$"), colnames(pca_summary))]*100,1),"%)")) + 
    theme_bw() + ggtitle(default.title) + my_theme_2 +
    theme(panel.grid.minor = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1) +
    scale_color_manual(values=pal)
    
  
  if(labels) {
    p <- p + ggrepel::geom_label_repel(aes(label = rownames(pca_plot), col=repel_col_by),
                              fontface = 'bold'
                              # , color = 'black'
                              , show.legend = F
                              , size=2,
                              box.padding = 0.35, point.padding = 0.5,
                              segment.color = 'grey50') 
  }
  
  if(scree_plot) {
    
    if(scree_plot_type=='explained_variance') {
      spdata      <- reshape2::melt(pca_summary[c(2:3),])
      spdata$Var2 <- as.numeric(gsub("PC", "",spdata$Var2))
      spdata      <- spdata[spdata$Var2<=max.pcs,]
      sp <- ggplot(spdata, aes(x=Var2, y=value, group=Var1, col=Var1)) + geom_point(size=1) + geom_line() +
        theme_bw() + my_theme_2 + scale_x_continuous(breaks = spdata$Var2) + 
        geom_hline(yintercept = 0.9, lwd=0.2, linetype="dashed", col = "darkred") +
        xlab("Principal Component") + ylab("Explained Variance") + ggtitle("PCA - Scree Plot") + 
        theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank()) +
        scale_y_continuous(breaks = seq(0,1,0.1),labels = scales::percent, limits = c(0,1)) + scale_color_manual(values=c("black","grey"))
    } else if(scree_plot_type=='standard_deviation') {
      spdata      <- reshape2::melt(pca_summary[1,])
      
      spdata$Var2 <- as.numeric(gsub("PC", "",rownames(spdata)))
      spdata      <- spdata[spdata$Var2<=max.pcs,]  
      
      sp <- ggplot(spdata, aes(x=Var2, y=value), col='black') + geom_point(size=1) + geom_line() +
        theme_bw() + my_theme_2 + scale_x_continuous(breaks = spdata$Var2) + 
        geom_hline(yintercept = 0.9, lwd=0.2, linetype="dashed", col = "darkred") +
        xlab("Principal Component") + ylab("Standard Deviation") + ggtitle("PCA - Scree Plot") + 
        theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank())
      
    }
    
    
    
    return(list("pca" = p, "screeplot" = sp))
  } else {
    return(p) 
  }
}

getPCA <- function(x) {
  
  rna <- t(x)
  rna <- rna[,colSums(rna)>1]
  RNA_pca <- prcomp(rna, scale. = T, center = T)
  pca_summary <- summary(RNA_pca)$importance
  
  return(list("RNA_pca"     = RNA_pca,
              "pca_summary" = pca_summary))
}


plot_loadings <- function(pca_data, components, assigned = NULL)
{
  if(is.list(pca_data)) {
    pca_load <- pca_data$RNA_pca$rotation
  } else {
    pca_load <- pca_data
  }
  
  toplot <- as.data.frame(pca_load[,components])
  if(!is.null(assigned)) {
    aidx <- match(assigned$variable, rownames(toplot))
    toplot$assigned <- assigned[aidx, 'princomp']
    nidx <- which(!toplot$assigned%in%components)
    toplot$assigned[nidx] <- "others"
    
    p0 <- ggplot(toplot, aes_string(x=components[1], y=components[2], col="assigned"))
  } else {
    p0 <- ggplot(toplot, aes_string(x=components[1], y=components[2]))
  }
  
  p <- p0 + geom_point(size=1, alpha=0.8) + 
    theme_bw() + my_theme + 
    geom_hline(yintercept = 0, lwd=0.2, linetype="dashed", col = "black") +
    geom_vline(xintercept = 0, lwd=0.2, linetype="dashed", col = "black") +
    ggtitle("PCA - Loading Plot") + scale_color_manual(values = c("#404040","#CC0000","orange")) +
    theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank()) 
  
  return(p)
}
assign_loadings <- function(pca_data, retained = NULL)
{
  if(is.list(pca_data)) {
    pca_load <- pca_data$RNA_pca$rotation
  } else {
    pca_load <- pca_data
  }
  
  if(!is.null(retained)) {
    if( is.numeric(retained) )  {
      if( length(retained)==1 )  retained <- 1:retained  
      retained <- paste0("PC",retained)
    }
    pca_load <- pca_load[,retained]
  }
  
  assigned <- apply(pca_load, 1, function(x) which.max(abs(x)))
  assigned <- paste0("PC",assigned)
  assigned <- data.frame("variable"    = rownames(pca_load)
                         , "princomp"  = assigned
                         , stringsAsFactors = F)
  return(assigned)
}
# Correlation Analysis ====
plotCorrelation <- function(x
                            , method = "pearson"
                            , myPalette=NULL
                            , ...) {
  
  mcor <- cor(x, method=method, ...)

  if(is.null(myPalette)) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))
  }
  hname <- paste0(method, " correlation")
  cHM <- ComplexHeatmap::Heatmap(mcor
                 , col  = myPalette(6)
                 , cell_fun = function(j, i, x, y, w, h, col) {
                         grid.text(round(mcor[i, j], digits = 2), x, y,gp = gpar(col='black', fontsize=6))
                       }
                 , name = hname
                 , row_names_gp = gpar(fontsize = 8)
                 , column_names_gp = gpar(fontsize = 8)
                 , heatmap_legend_param = list(title_position = "topcenter",
                                               # legend_width  = unit(4, "cm"),
                                               # legend_height = unit(0.5, "mm"),
                                               values_gp     = gpar(fontsize=8),
                                               legend_direction = "horizontal"))
  return(ComplexHeatmap::draw(cHM, heatmap_legend_side = "bottom"))
}

# Differential Expression Analysis ====
# still here for compatibility reason, 
# will be deprecated as soon as possible
calculateDiffExpr <- function(m
                              , group  = NULL
                              , filter = T
                              , filter.cpm.th = 1
                              , filter.sample.th = 2
                              , method = "exact"
                              , reference = NULL # Control group
                              , design = NULL # Custom design matrix
                              , cf = NULL # Testing coefficient (default: all columns in design matrix against control)
                              , contrast = NULL # Contrast matrix
                              , return.y  = F) {
  
  # Compute differential expression with EdgeR
  message("[*] Run EdgeR for Differential Expression Analysis")
  
  if(is.null(group)) {
    group <- as.factor(colnames(m))
    if(sum(table(levels(group))>1)==0) {
      stop(message("[!] Invalid experimental groups (< 2 replicates per condition)."))
    }
  } else {
    m <- m[,names(group)]
    if (!is.factor(group))   group <- as.factor(group)
    if (!is.null(reference)) group <- relevel(group, reference)
  }
  
  message(" -- Condition: ", paste0(levels(group), collapse = "-"))
  
  y <- edgeR::DGEList(counts=m, genes=rownames(m), group = group)
  y <- edgeR::calcNormFactors(y)
  
  # Clean environment
  rm(m)
  gc(verbose = F)
  
  if(filter) {
    message(" -- Filtering lowly expressed genes")
    message("    -- Threshods:")
    message("     * CPM = "     , filter.cpm.th)
    message("     * Samples = " , filter.sample.th)
    
    if(filter.sample.th>=1) {
      # Number of samples
      fs <- filter.sample.th
    } else {
      # Percentage of samples
      fs <- floor(ncol(m)*filter.sample.th)
    }
    
    keep <- rowSums(edgeR::cpm(y)>filter.cpm.th) >= fs
    y <- y[keep, , keep.lib.sizes=FALSE]
  }
  
  if(is.null(design)) {
    # Standard design matrix 
    # First column is control ( = reference in group)
    message(" -- Standard design matrix")
    design <- model.matrix(~group)
    colnames(design)[-1] <- paste0(levels(group)[-1], "vs", levels(group)[1])
  } else {
    # Customized design matrix
    message(" -- Custom model matrix")
  }
  
  if(is.null(cf) & is.null(contrast)) {
    # Test all columns in design matrix (ANOVA-like)
    # parametrize model wrt control ( = reference in group)
    cf  <- 2:ncol(design)
  }
  
  rownames(design) <- colnames(y)
  y <- edgeR::estimateDisp(y, design, robust=TRUE)
  
  message(" -- Testing differential expression, method: ", method)
  if(length(cf)>1) {
    message(" -- ANOVA-like for multiple group comparison") 
  }
  
  if(method=="exact") {
    # Exact test
    de <- edgeR::exactTest(y)
    de <- edgeR::topTags(de, n = Inf)
    
  } else if(method=="lrt") {
    # Likelihood-ratio test
    fit <- edgeR::glmFit(y, design)
    if(is.null(contrast)){
      # Standard comparison
      lrt <- edgeR::glmLRT(fit, coef=cf)
    } else {
      # GLM with contrasts
      lrt <- edgeR::glmLRT(fit, contrast = contrast)
    }
    de  <- edgeR::topTags(lrt, n = Inf)
    
  } else if(method=="qlf") {
    # Quasi-likelihood F-test
    fit <- edgeR::glmQLFit(y, design)
    if(is.null(contrast)){
      # Standard comparison
      qlf <- edgeR::glmQLFTest(fit, coef=cf)
    } else {
      # GLM with contrasts
      qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
    }
    de  <- edgeR::topTags(qlf, n = Inf)
    
  } else {
    stop(message("[-] Method not available. Please, provide a valid one ( exact / lrt / qlf )."))
  }
  
  if(return.y) {
    
    res <- list("y"  = y,
                "de" = de)
    return(res)
    
  } else {
    return(de)
  }
}

#' Differential Expression Analysis with edgeR
#'
#' Implementation of statistical models and tests 
#' for differential expresion analysis using edgeR
#' software package
#' @param y count matrix (or data.frame) or object of class DGEList
#' @param experimental_info data.frame containing informations about the experimental conditions for each sample
#' @param group conditions to test
#' @param method statistcal test to be performed (exact, lrt, qlf)
#' @return results of the test (data.frame)
#'
#' @examples
#' \dontrun{
#'     dea <- calculateDiffExprEdgeR(y, expreimental_info
#'     , group = experimental_info$group)
#' }
#' 
#'
#' @export
calculateDiffExprEdgeR <- function(y, experimental_info = NULL
                                   , gene_info = NULL
                                   , group  = NULL
                                   , reference = NULL # Control group
                                   , filter = T
                                   , filter.expr.th = 1
                                   , filter.sample.th = 2
                                   , normalize.using = "cpm"
                                   , tlen = NULL
                                   , method = "exact"
                                   , anovalike = F # Activate ANOVA-like test for any difference (alternative to cf,contrast)
                                   , design = NULL # Custom design matrix
                                   , cf = NULL # Testing coefficient (default: last column in design matrix)
                                   , contrast = NULL # Contrast matrix - character/numeric vector
                                   , return.y  = F) {
  
  # Compute differential expression with edgeR
  # Internals ---
  get_contrast <- function(contrast, design)
  {
    if(is.character(contrast) & length(contrast) == 3) {
      
      contrast <- limma::makeContrasts(paste0(contrast[1], contrast[-1], collapse = "-"), levels = design)
      
    } 
    return(contrast)
  }
  message("[*] Run edgeR for Differential Expression Analysis")
  
  if(is.matrix(y) | is.data.frame(y)) {
    y <- processRNAseqEdgeR(m = y
                            , experimental_info = experimental_info
                            , gene_info = gene_info
                            , group = group
                            , reference = reference
                            , filter = filter
                            , filter.expr.th = filter.expr.th
                            , filter.sample.th = filter.sample.th
                            , normalize.using = normalize.using
                            , tlen = tlen)
  }
  
  if(is.null(design)) {
    # Standard design matrix 
    # First column is control ( = reference in group)
    message(" -- Standard design matrix")
    design <- model.matrix(~group, data = y$samples)
    colnames(design)[-1] <- paste0(levels(group)[-1], "vs", levels(group)[1])
  } else {
    # Customized design matrix
    # Passed as string formula or model.matrix
    message(" -- Custom model matrix")
    if(is.character(design)) {
      design <- model.matrix(as.formula(design), data = y$samples)
    }
  }
  
  if(is.null(cf) & is.null(contrast)) {
    if(anovalike) {
      # Test all columns in design matrix (ANOVA-like)
      # parametrize model wrt control ( = reference in group)
      cf  <- 2:ncol(design)
    } else {
      # Default: last column in design matrix
      cf  <- ncol(design)
    }
  }
  
  rownames(design) <- colnames(y)
  y <- edgeR::estimateDisp(y, design, robust=TRUE)
  
  message(" -- Testing differential expression, method: ", method)
  if(length(cf)>1) {
    message(" -- ANOVA-like for multiple group comparison") 
  }
  
  if(method=="exact") {
    # Exact test
    de <- edgeR::exactTest(y)
    de <- edgeR::topTags(de, n = Inf)
    
  } else if(method=="lrt") {
    # Likelihood-ratio test
    fit <- edgeR::glmFit(y, design)
    if(is.null(contrast)){
      # Standard comparison
      lrt <- edgeR::glmLRT(fit, coef=cf)
    } else {
      # GLM with contrasts
      contrast <- get_contrast(contrast = contrast, design = design)
      lrt <- edgeR::glmLRT(fit, contrast = contrast)
    }
    de  <- edgeR::topTags(lrt, n = Inf)
    
  } else if(method=="qlf") {
    # Quasi-likelihood F-test
    fit <- edgeR::glmQLFit(y, design)
    if(is.null(contrast)){
      # Standard comparison
      qlf <- edgeR::glmQLFTest(fit, coef=cf)
    } else {
      # GLM with contrasts
      contrast <- get_contrast(contrast = contrast, design = design)
      qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
    }
    de  <- edgeR::topTags(qlf, n = Inf)
    
  } else {
    stop(message("[-] Method not available. Please, provide a valid one ( exact / lrt / qlf )."))
  }
  
  if(return.y) {
    res <- list("y"  = y,
                "de" = de)
    return(res)
    
  } else {
    return(de)
  }
}

saveXLSresEdgeR <- function(res, outfile, name, force=T) {
  # set java memory
  options(java.parameters = "-Xmx8000m")
  # java garbage collector
  # jgc <- function()
  # {
  #   .jcall("java/lang/System", method = "gc")
  # }
  
  message("[+] Saving results to file: ", outfile)
  
  if( is.list(res) & !is.data.frame(res)) {
    for( i in 1:length(res) ) {
      
      gc(verbose = F)
      # jgc()
      tox <- res[[i]]
      if(i==1 & file.exists(outfile) & force==F) {
        ap <- T
        warning(" -- appending to existing file.")
      } else {
        ap <- i>1
      }
      xlsx::write.xlsx2(tox, file = outfile, sheetName = names(res)[i], append = ap, row.names = F)
    }
  } else {
    # ap  <- F
    if(!is.data.frame(res)) res <- as.data.frame(res)
    tox <- res
    if(missing(name)) name <- "edgeR"
    if(file.exists(outfile) & force==F) {
      ap <- T
      warning(" -- appending to existing file.")
    } else {
      ap <- F
    }
    xlsx::write.xlsx2(tox, file = outfile, sheetName = name, append = ap, row.names = F)
  }
  
}

getDEres <- function(x, genes) {
  # Get DE results for specific genes
  if(!is.data.frame(x)) x <- x$table
  x[intersect(genes, rownames(x)),]  
}

getDEgenes <- function(x, fdrTh=0.1, fcTh=0.5, lcpmTh=NULL) {
  # Get DE genes satisfying thresholds
  message("[+] Get differentially expressed genes, thresholds:")
  message(" - logFC = ", fcTh)
  message(" - FDR = ", fdrTh)
  message(" - logCPM = ", lcpmTh)
  
  if(!is.data.frame(x)) x <- x$table
  
  if( !is.null(fdrTh) ) {
    x <- x[ x[,'FDR'] <= fdrTh , ]
  }
  
  if( !is.null(lcpmTh) ) {
    x <- x[ x[,'logCPM'] >= lcpmTh , ]
  }
  
  if(nrow(x) == 0) {
    stop(message("[-] No differentially expressed genes with provided FDR"))
  }
  
  if( !is.null(fcTh) ) {
    fcidx <- grep("logFC", colnames(x))
    
    if( length(fcidx)==1 ) {
      x <- x[ abs(x[,fcidx]) >= fcTh , ]
    } else {
      x <- x[ rowSums( abs(x[,fcidx]) >= fcTh) >= 1,]
    }
    
    if(nrow(x) == 0) {
      stop(message("[-] No differentially exxpressed genes with provided logFC"))
    }
  }
  
  return(x)
  
}

getDEgsigned <- function(deg, signed=(-1))
{
  fcidx <- grep("logFC", colnames(deg))
  x     <- deg[,fcidx]
  
  max.idx <- apply(x, 1, function(m) which.max(abs(m)))
  
  y <- vector(mode = 'numeric', length = nrow(x))
  names(y) <- names(max.idx)
  
  for(i in names(max.idx)) {
    y[i] <- x[i, max.idx[i]]
  }
  
  gene.idx <- names(y[which(sign(y)==signed)])
  
  return(x[gene.idx,])
}


plotRNAVolcanos <- function(de, lfcTh=1, pvTh=0.05, top=5, gtitle = NULL)
{
  fcIdx <- grep("logFC|log2FoldChange", colnames(de), value = T)
  pvIdx <- grep("adj.P.Val|FDR", colnames(de), value = T)
  tmp   <- de[,c(fcIdx,pvIdx)]
  colnames(tmp) <- c("lfc","padj")
  tmp$status <- "none"
  tmp$status[which(tmp$lfc>=lfcTh & tmp$padj<=pvTh)] <- "Up-regulated"
  tmp$status[which(tmp$lfc<=(-lfcTh) & tmp$padj<=pvTh)] <- "Down-regulated"
  tmp$status <- factor(tmp$status)
  
  topUP <- tmp[tmp[,'status']=="Up-regulated",]
  nup <- nrow(topUP)
  topDW <- tmp[tmp[,'status']=="Down-regulated",]
  ndw <- nrow(topDW)
  
  if(!is.null(top)) {
    topUP <- rownames(topUP[1:top,])
    topDW <- rownames(topDW[1:top,])
    
    tmp$lab <- rownames(tmp)
    tmp$lab[-which(rownames(tmp)%in%c(topDW, topUP))] <- ""
  } else {
    tmp$lab <- ""
  }
  
  p <- ggplot(tmp, aes(x=lfc, y=-log10(padj),col=status, label=lab)) + geom_point(size=1, alpha=0.8) + theme_bw() +
    geom_hline(yintercept = -log10(pvTh), linetype = 'dashed', lwd = 0.25) + 
    geom_vline(xintercept = lfcTh, linetype = 'dashed', lwd = 0.25) +
    geom_vline(xintercept = -lfcTh, linetype = 'dashed', lwd = 0.25) +
    xlab("logFC") + ylab("adjusted P-value") +
    theme_bw() + my_theme + ggtitle(gtitle) +
    theme(plot.title = element_text(size=10, face = "bold", hjust = 0.5)) +
    scale_color_manual(values = c('#004C99','#404040','#CC0000')) +
    ggrepel::geom_text_repel(sshow.legend = F) +
    geom_label(
      data    = subset(tmp, status=="Up-regulated"),
      mapping = aes(x = max(tmp$lfc)-1.5, y = max(-log10(tmp$padj))-30, label = nup, col = status),
      size = 3
    ) +
    geom_label(
      data    = subset(tmp, status=="Down-regulated"),
      mapping = aes(x = min(tmp$lfc)+1.5, y = max(-log10(tmp$padj))-20, label = ndw, col = status),
      size = 3
    ) 
  
  return(p)
}

build_model_matrix <- function(group, ref, formula_string)
{
  # build model --
  if(is.character(group) & !missing(ref)) {
    # parametrize with "ref" as baseline 
    group <- factor(group)
    group <- relevel(group, ref = ref)
    mod   <- model.matrix(~group)
    colnames(mod)[-1] <- paste0(levels(group)[-1], "vs", levels(group)[1])
  } else {
    if(missing(formula_string)) {
      stop(message("[!] Please, provide parametrization formula (as charachter)"))
    }
    fm  <- as.formula(formula_string)
    mod <- model.matrix(fm, data = group)
  }
  
  return(mod)
}

prepare_pairwise_contrasts <- function(design, group, exclude = NULL)
{
  if(missing(design)) {
    group  <- as.factor(group)  
    design <- model.matrix(~group)
  }
  
  pairwise  <- apply(combn(rev(colnames(design)), 2),2, paste, collapse = "-")
  pairwise  <- gsub("\\-\\(Intercept\\)","", pairwise)
  if(!is.null(exclude)) {
    pairwise <- pairwise[grep(exclude, pairwise, invert = T)]
  }
  contrasts <- limma::makeContrasts(contrasts = pairwise, levels = design)
  
  return(list("contrasts" = contrasts, "design" = design))
}

fix_de_contrast_names <- function(de, ref)
{
  idx_nm <- grep("-", names(de))
  tmp_nm <- names(de)[idx_nm]
  names(de)[idx_nm] <- gsub("-","vs",gsub(paste0("vs",ref),"",tmp_nm))
  
  return(de)
}

analyze_de_contrasts <- function(counts, design, test_contrasts
                                 , method  = "lrt"
                                 , fdrTh   = 0.05
                                 , fcTh    = 1
                                 , lcpmTh  = 0.5
                                 , fixnm_ref = NULL
                                 , ...)
{
  
  de <- lapply(colnames(test_contrasts), function(cnt) {
    message("Contrast: ", cnt)
    x <- calculateDiffExpr(m          = counts
                           , design   = design
                           , contrast = test_contrasts[,cnt]
                           , method   = method
                           , ...)
    y <- getDEgenes(x, fdrTh = fdrTh, fcTh = fcTh, lcpmTh = lcpmTh)
    res <- list("table" = x$table, "sig" = y)
    return(res)
  })
  
  names(de) <- colnames(test_contrasts)
  
  if(!is.null(fixnm_ref)) de <- fix_de_contrast_names(de, ref = fixnm_ref)
  
  return(de)
}

create_fc_matrix_contrasts <- function(de_contrasts)
{
  gidx <- lapply(de_contrasts, function(i) i$sig$genes)
  gidx <- unique(unlist(gidx))
  
  fc <- lapply(de_contrasts, function(i) i$table[gidx,"logFC"])
  fc <- do.call(cbind, fc)
  rownames(fc) <- gidx
  
  return(fc)
}

# Wrapper DE analysis
analyzeDE <- function(s
                      , fdrTh  = 0.05
                      , fcTh   = 0.5
                      , lcpmTh = 0
                      , multigroup_method = "qlf"
                      , ...)
{
  method <- NA
  deAnalysis    <- lapply(s, function(x) {
    if(length(unique(x))>2) {
      method <<- multigroup_method
    } else {
      method <<- "exact"
    }
    tryCatch( calculateDiffExpr(m=counts, group = x, method = method, ...)
              , error = function(e) {
                message(paste0("[-] Error: ",e, ".\n"))
                return(NA)
              }  )
  } )
  
  deAnalysis.Sig <- lapply(deAnalysis, function(x) {
    tryCatch(getDEgenes(x, fdrTh = fdrTh, fcTh = fcTh, lcpmTh = lcpmTh)
             , error   = function(e) return(NA))
  })
  
  nDeGenes  <- lapply(deAnalysis.Sig, nrow)
  
  
  if(method == "qlf") {
    up <- lapply(deAnalysis.Sig, function(x) {
      tryCatch(getDEgsigned(x, signed=1)
               , error   = function(e) return(NA))
    })
    
    dw <- lapply(deAnalysis.Sig, function(x) {
      tryCatch(getDEgsigned(x, signed=(-1))
               , error   = function(e) return(NA))
    })
  } else if(method == "exact"){
    up <- lapply(deAnalysis.Sig, function(x) {
      tryCatch(x[sign(x[,'logFC'])==1,]
               , error   = function(e) return(NA))
    })
    
    dw <- lapply(deAnalysis.Sig, function(x) {
      tryCatch(x[sign(x[,'logFC'])==(-1),]
               , error   = function(e) return(NA))
    })
  }
  
  
  message("[+] Plotting ...")
  
  hm <- lapply(names(s), function(i) {
    tryCatch({deg <- rownames(deAnalysis.Sig[[i]])
    gidx <- intersect(rownames(tpm), deg)
    
    set.seed(s33d)
    m <- tpm[gidx, names(s[[i]])]
    h <- get_heatmap3(m
                      , show_row_names = F
                      , clustering_distance_rows = "euclidean"
                      , clustering_distance_columns = "euclidean"
                      , clustering_method_rows = "complete"
                      , clustering_method_columns = "complete"
                      , show_row_dend = F
                      , retHm = T)
    return(h)}, error   = function(e) return(NA))
  })
  names(hm) <- names(s)
  
  vp <- lapply(deAnalysis, function(i) plotRNAVolcanos(de = i$table, lfcTh = fcTh, pvTh = fdrTh))
  
  return(list(    "deAnalysis"     = deAnalysis
                  , "deAnalysis.Sig" = deAnalysis.Sig
                  , "nDeGenes"       = nDeGenes
                  , "up"             = up
                  , "dw"             = dw
                  , "heatmap"        = hm
                  , "vp"             = vp))
}

saveDE <- function(dea, FIGDIR, RESDIR)
{
  message("[+] Saving plots, directory: ", FIGDIR)
  message(" -- volcanos ")
  lapply(names(dea$vp), function(i) {
    outfig <- paste0(FIGDIR, "/vp.analysis_",i
                     ,".fcTh_",fcTh
                     ,".fdrTh_",fdrTh
                     ,".pdf")
    pdf(file = outfig, paper = 'a4', width = unit(6,'cm'), height = unit(6,'cm'))
    print(dea$vp[[i]])
    dev.off()
  }
  )
  message(" -- heatmaps ")
  lapply(names(dea$heatmap), function(i) {
    outfig <- paste0(FIGDIR, "/hm.analysis_",i
                     ,".fcTh_",fcTh
                     ,".fdrTh_",fdrTh
                     ,".pdf")
    set.seed(s33d)
    pdf(file=outfig, paper = 'a4', width = unit(3,'cm'), height = unit(6,'cm'))
    set.seed(s33d); print(dea$heatmap[[i]])
    dev.off()
  }
  )
  
  message(" -- tables ")
  if(!exists("analysis")) analysis <- "analysis"
  outfile <- paste0(RESDIR,"/DEA.edgeR.",analysis
                    ,".fcTh_",fcTh
                    ,".fdrTh_",fdrTh
                    ,".xlsx")
  
  saveXLSresEdgeR(dea$deAnalysis.Sig, outfile = outfile)
}

compute_fold_change <- function(x, y, pseudocount, logscale = F) 
{
  if(logscale) {
    fc <- x-y
  } else {
    fc <- log2((x+pseudocount)/(y+pseudocount))  
  }
  return(fc)
} 

analyze_relative_expression <- function(genes, dea, time_0 = "T0h", facet_formula = NULL)
{
  get_expression_fc <- function(dea, genes)
  {
    summarize_fc <- function(x, direction)
    {
      nx <- nrow(x)
      x <- as.data.frame(t(apply(x, 2, data_summary)))
      x$time <- rownames(x)
      x$direction <- direction
      x$n <- nx
      
      return(x)
    }
    
    if(!missing(genes)) {
      sidx    <- intersect(genes, dea$sig$genes)
      dea$sig <- dea$sig[sidx,]
      dea$table <- dea$table[genes,]
    }
    up <- getDEgsigned(dea$sig, signed = 1)
    gup <- rownames(up)
    up <- summarize_fc(x = up, direction = "up")
    dw <- getDEgsigned(dea$sig, signed = (-1))
    gdw <- rownames(dw)
    dw <- summarize_fc(x = dw, direction = "dw")
    
    cgenes <- setdiff(dea$table$genes, dea$sig$genes)
    const  <- dea$table[cgenes, grep("FC", colnames(dea$table))]
    gconst <- rownames(const)
    const <- summarize_fc(x = const, direction = "const")
    
    toplot <- rbind.data.frame(up, dw, const)  
    toplot <- plyr::unrowname(toplot)
    toplot$time <- gsub("logFC.","", toplot$time)
    toplot$time <- sapply(strsplit(toplot$time, "_"),"[[",1)
    
    t0 <- data.frame("y" = rep(0,3)
                     , "ymin" = rep(0,3)
                     , "ymax" = rep(0,3)
                     , "time" = rep(time_0,3)
                     , "direction" = c("up","dw","const")
                     , "n" = c(unique(up$n),unique(dw$n),unique(const$n)))
    
    toplot <- rbind.data.frame(toplot, t0)  
    
    return(list("toplot"=toplot, "gup" = gup, "gdw" = gdw, "gconst" = gconst))
  }
  wrap_get_efc <- function(dea, genes)
  {
    efc <- lapply(dea, function(i) 
    {
      efc <- get_expression_fc(i, genes)
      
    })
    
    reg <- lapply(efc, function(i) {
      i$toplot <- NULL
      i
    } )
    
    efc <- lapply(efc, "[[", "toplot")
    efc <- do.call(rbind, efc)
    efc$cond <- gsub("\\.+\\d+", "",rownames(efc))
    
    return(list("efc"=efc,"reg"=reg))
  }
  plot_expression_fc <- function(efc, facet_formula = NULL)
  {
    if(is.list(efc) & !is.data.frame(efc)) {
      efc <- do.call(rbind, efc) 
      efc$group <- sapply(strsplit(rownames(efc), "[.]"),"[[",1)
      if(is.null(facet_formula)) {
        p0 <- ggplot(efc, aes(x=time, y=y)) + facet_grid(group~cond) 
      } else {
        p0 <- ggplot(efc, aes(x=time, y=y)) + facet_wrap(as.formula(facet_formula), ncol = 3) 
      }
    } else {
      p0 <- ggplot(efc, aes(x=time, y=y)) + facet_grid(~cond) 
    }
    
    p <- p0 +
      theme_bw() + my_theme + ylab("average logFC") +
      geom_line(aes(group = direction, col=direction), size = 0.5) +
      geom_label(
        data    = subset(efc, direction=="up"),
        mapping = aes(x = 1, y = 3, label = n, col = direction),
        size = 3
      ) +
      geom_label(
        data    = subset(efc, direction=="dw"),
        mapping = aes(x = 1, y = -2, label = n, col = direction),
        size = 3
      ) +
      geom_ribbon(aes(ymin=ymin, ymax=ymax, group=direction, fill=direction), alpha = 0.2) +
      scale_fill_manual(values = c("lightgrey","blue","red")) +
      scale_color_manual(values = c("lightgrey","blue","red")) +
      theme(panel.grid = element_blank()
            , axis.title.x = element_blank()
            , axis.text.x = element_text(angle = 90, hjust = 1)
            , strip.background = element_blank()
            , legend.position =  "none"
            # , strip.text = element_text(face="bold", size = 8)
      ) 
    
    return(p)
  }
  
  if(is.list(genes)) {
    tmp <- lapply(genes, wrap_get_efc, dea = dea)
    efc <- lapply(tmp, "[[", "efc")
    reg <- lapply(tmp, "[[", "reg")
    rm(tmp)
    efc <- list("efc"=efc,"reg"=reg)
  } else {
    efc <- wrap_get_efc(dea = dea, genes = genes)
  }
  
  p <- plot_expression_fc(efc$efc, facet_formula = facet_formula)
  
  return(list("p" = p, "data" = efc))
} 

# DESeq2 ---
calculateDiffExprDESeq2 <- function(counts
                                    , info_analysis
                                    , design_formula 
                                    , contrast_tests = NULL
                                    , fcth           = NULL
                                    , pvth           = NULL
                                    , outfile        = NULL) 
{
  message("[*] Run DESeq2 for Differential Expression Analysis")
  m           <- counts[,rownames(info_analysis)]
  m           <- apply(m, 2, as.integer)
  rownames(m) <- rownames(counts)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m,
                                colData   = info_analysis,
                                design    = as.formula(design_formula))
  dds <- DESeq2::DESeq(dds)
  
  if(!is.null(contrast_tests)) {
    res <- list()
    for(i in contrast_tests) {
      k <- paste0(i[1],"_",i[2],"_vs_",i[3])
      if(!is.null(fcth) & !is.null(pvth)) {
        res[[k]] <- list()
        res[[k]]$result <- DESeq2::results(dds, contrast = i)
        res[[k]]$sig <- subset(DESeq2::results(dds, contrast = i), abs(log2FoldChange)>=fcth & padj <= pvth)
        res[[k]]$sig <- cbind.data.frame("genes"=rownames(res[[k]]$sig),res[[k]]$sig)
        res[[k]]$sig <- res[[k]]$sig[order(res[[k]]$sig$padj, decreasing = F),]
      } else {
        res[[k]] <- DESeq2::results(dds, contrast = i)
      }
    }
    if(!is.null(outfile)) {
      saveXLSresEdgeR(lapply(res,"[[","sig"), outfile = outfile, force = T)
    }
    
    return(list("dds"=dds,"res"=res))
  } else {
    return(dds)
  }
}

# Surrogate Variable Analysis ====
cleaningP <- function(y, mod, svaobj, P=ncol(mod))
{
  # adapted from Jaffe et al. (BMC Bioinformatics 2015):
  # Practical impacts of genomic data “cleaning” on 
  # biological discovery using surrogate variable analysis
  
  X      <- cbind(mod,svaobj$sv)
  Hat    <- solve(t(X)%*%X)%*%t(X)
  beta   <- (Hat%*%t(y))
  
  cleany <- y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  
  return(cleany)
}
estimate_surrogate_vars <- function(y
                                    , mod = NULL # Custom model
                                    , clean = F
                                    , norm.var = "CPM") 
{
  if(class(y)!="DGEList") {
    stop(message("[!] Please, provide DGEList object."))
  }

  # build model --
  if(is.null(mod)) {
    group <- relevel(y$samples$group, ref = y$ref)
    mod   <- model.matrix(~group)
    colnames(mod)[-1] <- paste0(levels(group)[-1], "vs", levels(group)[1])
  }
  
  # estimate surrogate variables --
  svseq <- sva::svaseq(y[[norm.var]], mod, mod[,1])
  colnames(svseq$sv) <- paste0("sv", 1:ncol(svseq$sv)) 
  mod.sv <- cbind(mod, svseq$sv)
  
  res <- list("mod.sv" = mod.sv, "svseq" = svseq)
  
  if(clean) {
    # regress out surrogate variables from data --
    y$cidx <- cleaningP(y = y[[norm.var]], mod = mod, svaobj = svseq)
    names(y)[grep('cidx', names(y))] <- paste0('clean', norm.var) 
    
    res <- c(list("y" = y), res)
  }
  
  return(res)
}

# Clustering ====
# K-means ---
findClustSSE <- function(scaledata, nKM=20, ret=F)
{
  wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
  for (i in 2:nKM) wss[i] <- sum(kmeans(scaledata,
                                       centers=i)$withinss)
  wss.best <- as.numeric(which.max(wss))
  plot(1:nKM, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  abline(v = wss.best, lty = 2)
  
  if(ret) return(wss.best)
}

findClustASW <- function(scaledata, nKM=20, nstart = 100, ret=F)
{
  sil <- rep(0, nKM)
  #repeat k-means for 1:n and extract silhouette:
  for(i in 2:nKM){
    k1to20 <- kmeans(scaledata, centers = i, nstart = nstart, iter.max = 20)
    ss <- silhouette(k1to20$cluster, dist(scaledata))
    sil[i] <- mean(ss[, 3])
  }
  
  # Plot the  average silhouette width
  plot(1:nKM, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
  abline(v = which.max(sil), lty = 2)
  
  sil.best <- as.numeric(which.max(sil))
  cat("Average silhouette width optimal number of clusters:", sil.best, "\n")
  
  if(ret) return(sil.best)
}

findClustCal <- function(scaledata, nKM=20, ret=F)
{
  # Calinski-Harabasz index
  fit <- vegan::cascadeKM(scaledata, 1, nKM, iter = 100)
  plot(fit, sortg = TRUE, grpmts.plot = TRUE)
  calinski.best <- as.numeric(which.max(fit$results[2,]))
  cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
  
  if(ret) return(calinski.best)
}

setKMClusters <- function(m, method="silhouette", ...)
{
  scaledata <- t(scale(t(m)))
  
  if( method=="sse" ) {
    return(findClustSSE(scaledata, ...))
  }
  
  if( method=="silhouette") {
    return(findClustASW(scaledata, ...))
  }
  
  if( method=="calinski") {
    return(findClustCal(scaledata, ...))
  } else {
    stop(message('[!] Invalid method, please provide one of sse/silhouette/calinksi'))
  }
}

# partitioning aroung medoids (PAM) ---
setPAMCnumber <- function(m, metric = "euclidean", scaledata = T, nc = 12, ret = T)
{
  if(scaledata) m <- t(scale(t(m)))
  pams <- lapply(2:nc, function(k) cluster::pam(m, k = k, metric = metric))
  names(pams) <- 2:nc
  sw <- unlist(lapply(pams, function(i) i$silinfo$avg.width))
  oc <- as.numeric(names(pams)[which.max(sw)])
  message(" - PAM optimal number of clusters = ",  oc)
  
  if(ret) return(oc)
}

getPAM <- function(m, k, metric = "euclidean", scaledata = T, ...) {
  if(scaledata) m <- t(scale(t(m)))
  list('cluster' = cluster::pam(m,k, metric = metric, ...))
}

# Hierarchical clustering ---
get_hclust <- function(m
                       , distance = 'euclidean'
                       , method   = 'complete'
                       , cut_tree = 'static'
                       , cut_tree_h = NULL
                       , cut_tree_k = NULL
                       , return_d   = F
                       , return_hr  = T
                       , minClusterSize = 50
                       , ...)
{
  
  message("[*] Hierarchical clustering")
  message(" -- distance: ", distance)
  message(" -- method: ", method)
  
  m  <- t(scale(t(m), center = T, scale = T))
  
  res <- vector(mode = 'list')
  
  if( any(grepl("pearson|spearman|kendall", distance)) ) {
    c <- cor(t(y), method = distance) 
    d <- as.dist(1-c)
  } else {
    d <- dist(m, method = distance)
  }
  
  if(return_d) {
    res$d <- d
  }
   
  hr <- hclust(d, method = method, members=NULL)
  
  if(return_hr) {
    res$hr <- hr
  }
  
  if(cut_tree == 'static') {
    
    if(is.null(cut_tree_h) || is.null(cut_tree_k)) {
      cut_tree_h <- quantile(hr$height, probs = seq(0,1,0.01))['99%']  
    } 
    message(" -- Static tree cutting")
    message(" -- heigth = ", cut_tree_h)
    cl <- cutree(hr, h=cut_tree_h, cut_tree_k)
    
    res$cl <- cl
    rm(cl)
      
  } else if(cut_tree == 'dynamic') {
    
    message(" -- Dynamic tree cutting")
    
    cl <-  cutreeDynamic(dendro = hr
                         , distM = as.matrix(d)
                         , minClusterSize = minClusterSize
                         , deepSplit = 0)
    res$cl <- cl
    rm(cl)
    
  } else if(!return_d & !return_hr) {
    res$hr <- hr
  } 
  
  return(res)
}
# Heatmap ====
get_heatmap3 <- function(m
                         , annotDF  = NULL
                         , annotCol = NULL
                         , fig_out  = NULL
                         , retHm    = F
                         , bm = T
                         , myPalette = NULL
                         , myZscale = NULL
                         , myLegend = NULL
                         , scale = T
                         , ...){
  
  base_mean <- rowMeans(m)
  if(scale) {
    m_scaled <- t(apply(m, 1, scale))
    colnames(m_scaled) <- colnames(m)
  } else {
    m_scaled <- m
  }
  
  bPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))
  
  if(is.null(myPalette)) {
    myPalette <- c("blue","black","red")
  }
  
  if(is.null(myZscale)) {
    myZscale <- c(-2, 0, 2)
  }
  ramp <- circlize::colorRamp2(myZscale, myPalette)
  
  
  if (!is.null(annotDF)) {
    if (!is.null(annotCol)) {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF, 
                                     col = annotCol, 
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8)))
    } else {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF, 
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8)))
    }
  } else {
    ha_column <- new("HeatmapAnnotation")
  }
  
  if(is.null(myLegend)) myLegend <- "TPM" 
  
  hm <- ComplexHeatmap::Heatmap(m_scaled, col = ramp,
                # show_row_dend = T,
                # row_names_side = "left",
                # row_names_gp = gpar(fontsize=8),
                column_names_gp = gpar(fontsize=8),
                column_title_gp = gpar(fontsize=10, fontface="bold"),
                heatmap_legend_param = list(title = paste0("Z-score (",myLegend,")"),
                                            title_gp = gpar(fontsize=8),
                                            title_position = "topcenter",
                                            legend_width  = unit(3, "cm"),
                                            legend_height = unit(0.5, "mm"),
                                            values_gp     = gpar(fontsize=8),
                                            # legend_direction = "vertical"
                                            legend_direction = "horizontal"
                                            )
                , top_annotation = ha_column
                , top_annotation_height = unit(4, "mm")
                , ...)
  
  if(bm) {
    bmscale <- summary(base_mean)
    bmramp <- circlize::colorRamp2(c(bmscale[1],bmscale[3],bmscale[5]), bPalette(3))
    bmh <- ComplexHeatmap::Heatmap(base_mean 
                   # , name = "Mean Expression"
                   , column_names_gp = gpar(fontsize=8)
                   , show_row_names = FALSE
                   , width = unit(3, "mm") 
                   , col = bmramp
                   , heatmap_legend_param = list(title = paste0("Average ",myLegend),title_gp = gpar(fontsize=8)))
    
    hmOut <- hm + bmh
  } else {
    hmOut <- hm 
  }
  
  if(!is.null(fig_out)){
    pdf(file = fig_out, useDingbats = F, h=8, w=3, paper = "a4")
    ComplexHeatmap::draw(hmOut, heatmap_legend_side = "right")
    dev.off()
  } else{
    ComplexHeatmap::draw(hmOut, heatmap_legend_side = "right")
  }
  
  if(retHm) return(hmOut)
}

get_clusters <- function(m, hm){
  # Retrieve clusters from K-means
  clusters <- lapply(row_order(hm), 
                     function(x){
                       rownames(m[x,])
                     }
  )
  names(clusters) <- paste0('cluster_', seq_along(clusters))
  return(clusters)
}


plot_cluster_expression <- function(m, cl, pal)
{
  toplot <- lapply(cl, function(x) y <- m[x,])
  toplot <- reshape2::melt(toplot)
  colnames(toplot) <- c("gene","sample","tpm","cluster")
  toplot <- plyr::ddply(toplot, .(gene, sample, cluster)
                  , summarize
                  , av_tpm = mean(tpm)
                  , log2_av_tpm = log2(mean(tpm)))
  
  toplot$state <- gsub("/.*","",as.character(toplot$sample))
  toplot$cond  <- gsub(".*/| #.*","",as.character(toplot$sample))
  
  # toplot$sample <- factor(toplot$sample, levels = c("mES/WT #1","mES/WT #2","mES/3BKO #1","mES/3BKO #2"
  #                                                   ,"EpiSC/WT #1","EpiSC/WT #2","EpiSC/3BKO #1","EpiSC/3BKO #2") )
  toplot$cluster <- gsub("_"," ", toplot$cluster)
  p <- ggplot(toplot, aes(x=sample, y=log2_av_tpm, fill = cond)) + 
    geom_boxplot(notch = T, outlier.colour = "grey", outlier.size = 0.2, outlier.alpha = 0.5) +
    facet_wrap(~cluster, ncol = 1) + theme_bw() + my_theme + ylab("average log2[RPKM]") +
    # coord_cartesian(ylim = quantile(toplot$log2_av_tpm, c(0.1, 0.95))) +
    scale_fill_manual(values = pal) + 
    theme(panel.grid = element_blank()
          , axis.title.x = element_blank()
          , axis.text.x = element_text(angle = 90, hjust = 1)
          , strip.background = element_blank())
 
}

plot_hm_fc <- function(res_df
                       , myPalette = NULL
                       , myFCscale = NULL
                       , myTitle = NULL
                       , ...)
{
  
  
  if(any(grepl("FDR", colnames(res_df)))) {
    idx <- grep("logFC", colnames(res_df))
    fc <- as.matrix(res_df[,idx])
    colnames(fc) <- gsub("logFC.group|NT_|Act_","", colnames(fc))
  } else {
    fc <- res_df
  }
  if(is.null(myPalette)) {
    myPalette <- rev(colorRampPalette(RColorBrewer::brewer.pal(6, "RdBu"))(3))
  }
  
  if(is.null(myFCscale)) {
    myFCscale <- c(-1, 0, 1)
  }
  ramp <- circlize::colorRamp2(myFCscale, myPalette)
  
  if(is.null(myTitle)) {
    myTitle <- "log2[fold-change]"
  }
  hmfc <-  ComplexHeatmap::Heatmap(fc, col = ramp,
                   # show_row_dend = T,
                   # row_names_side = "left",
                   # row_names_gp = gpar(fontsize=8),
                   column_names_gp = gpar(fontsize=8),
                   column_title_gp = gpar(fontsize=10, fontface="bold"),
                   heatmap_legend_param = list(title = myTitle,
                                               title_gp = gpar(fontsize=8),
                                               title_position = "topcenter",
                                               legend_width  = unit(3, "cm"),
                                               legend_height = unit(0.5, "mm"),
                                               values_gp     = gpar(fontsize=8),
                                               # legend_direction = "vertical"
                                               legend_direction = "horizontal")
                   , ...)
  
  return(hmfc)
}
# Gene Ontology ====
getGO <- function(x, reg="UP", ORGANISM = "mmusculus", gl_input=T){
  
  # Perform GO term enrichment using gProfileR
  if(gl_input) {
    gene_list <- x
  } else {
    if ( reg!="both") {
      gene_list <- subset(x, regulation==reg)$symbol 
    } else {
      gene_list <- x$symbol
    }
  }
  
  gProfileR::gprofiler( gene_list,
             organism = ORGANISM,
             # custom_bg = background,
             max_p_value = 0.05,
             ordered_query =T,
             significant = T,
             correction_method='gSCS',
             min_set_size = 10,
             # evcodes=T,
             # exclude_iea = T,
             hier_filtering ='moderate',
             # src_filter = collections
             # png_fn = "Prostate/Figures/Goprofiler.AS.png"
             # include_graph=T
             src_filter = c("GO",'KEGG',"REAC","CORUM","OMIM",  "HP")
  )
  
}

processGO <- function(x, pvTh=0.05, cut=F){
  
  x <- subset(x, p.value<=0.05)
  x <- x[,c("term.name", "p.value")]
  if( any(duplicated(x$term.name))) x <- x[-which(duplicated(x$term.name)),]
  x$term.name <- factor(x$term.name, 
                        levels = x$term.name[order(x$p.value, decreasing = F)])
  if(cut){
    if(nrow(x)>25) x <- x[1:25,]
  }
  x <- x[with(x, order(p.value, decreasing = F)),]
  return(x)
}


plotGObars <- function(go, tool='clusterProfiler', ntop=NULL, ...)
  {
  
  if(!is.null(ntop)) {
    go <- go[1:ntop,]
  }
  
  if( tool == 'gprofiler') {
    ggplot(go, aes(x=term.name, y=-log10(p.value))) + geom_col(position = "dodge") +
      coord_flip() + theme_bw() + my_theme + theme(panel.grid = element_blank())
  } else if( tool == 'clusterProfiler') {
    go$Description <- factor(go$Description, 
                          levels = go$Description[order(go$qvalue, decreasing = T)])
    ggplot(go, aes(x=Description, y=-log10(qvalue))) + geom_col(position = "dodge", ...) +
      coord_flip() + theme_bw() + my_theme + theme(panel.grid = element_blank())
  }
}

getGO_v2 <- function(geneList
                     , species      = 'mm'
                     , custom       = F
                     , qvalueCutoff = 0.05
                     , ont          = "BP"
                     , simplify     = T
                     , maxGSSize    = 500
                     , ...)
{
  if(species=='mm') {
    require(org.Mm.eg.db)
    db <- org.Mm.eg.db
    
  } else if(species=='hg') {
    require(org.Hs.eg.db)
    db <- org.Hs.eg.db
    
  } 
  
  if(is.null(ont)) ont <- "BP"
  
  if(!custom){
    ego <- clusterProfiler::enrichGO(   gene         = geneList
                       , OrgDb         = db
                       , keyType       = 'SYMBOL'
                       , ont           = ont
                       , pAdjustMethod = "BH"
                       , pvalueCutoff  = 0.05
                       , qvalueCutoff  = qvalueCutoff
                       , maxGSSize     = maxGSSize)
  } else {
    ego <- clusterProfiler::enrichGO(  gene         = geneList
                      , OrgDb        = get(db)
                      , ... )
  }
  
  if(simplify) {
    ego <- simplify(ego)
  }
  
  return(ego)
}

get_kegg_pathway_enrichment <- function(geneList
                     , species      = 'mmu' # 'hsa'
                     , pvalueCutoff = 0.05
                     , method = 'enrichment' # 'gsea'
                     , ...)
{
  if(method=='enrichment') {
    kk <- clusterProfiler::enrichKEGG(  gene         = geneList
                     , organism     = species
                     , pvalueCutoff = pvalueCutoff
                     , ...)
    
  } else if(method=='gsea') {
    
    kk <- clusterProfiler::gseKEGG(  geneList     = geneList
                   , organism     = species
                   , nPerm        = 1000
                   , pvalueCutoff = pvalueCutoff
                   , verbose      = FALSE
                   , ...)
  }
  
  return(kk)
}

get_top_enrichment <- function(enrichment, top=10, ranking = "p.adjust")
{
  get_top <- function(x) 
  {
    if(grepl("DataFrame|data.frame",class(x))) {
      if(any(grepl("term.name",colnames(x)))) {
        colnames(x) <- c("Description","p.adjust")
        x$Count <- NA
      } else {
        x <- x[,c(1,4:5)] 
        colnames(x) <- c('Description','p.adjust','Count')
      }
      
    } else {
      if(grepl("enrichResult",class(x))) x <- x@result
      if(ranking=="p.adjust") {
        x <- x[with(x, order(p.adjust, decreasing = F)),]
      } else if(ranking=="Count"){
        x <- x[with(x, order(Count, decreasing = T)),]
      }
    }
    if(nrow(x)<top) top <- nrow(x)
    x[1:top,c('Description','p.adjust','Count')]
  }
  
  if(is.list(enrichment)) {
    y <- lapply(enrichment, get_top)
  } else{
    y <- get_top(x = enrichment)
  }
  return(y)
}


plot_top_enrichment <- function(top_enrichment, transpose = F, ...) 
{
  get_matrix <- function(x)
  {
    m <- as.matrix(x[,2])
    rownames(m) <- x[,1]
    return(m)
  }
  
  def.pal <- function(pal)
  {
    if(missing(pal)) {
      message(" -- using default palette (Reds).")
      pal  <- c('white', RColorBrewer::brewer.pal(9, "Reds"))
    } else {
      message(" -- using custom defined palette.")
    }
    return(pal)
  }
  
  plot.heatmap <- function(m
                           # , log_scale=F
                           , term_font_size = 8
                           ,...)
  {
    pal <- def.pal(...)
    ramp <- circlize::colorRamp2(seq(0,length(pal)-1,1), pal)
    # if(log_scale) {
    #   ramp <- colorRamp2(seq(-1,1,0.6), pal)
    #   idx <- m>0
    #   m[idx] <- log2(m[idx])
    #   }
    hm <- ComplexHeatmap::Heatmap(m
                  , show_column_names = T
                  , show_row_names    = T
                  , col = ramp
                  , cluster_rows = F
                  , cluster_columns = F
                  , heatmap_legend_param = list(title = "-log10[adjusted-p]",
                                                title_gp = gpar(fontsize=8),
                                                title_position = "topcenter",
                                                legend_width  = unit(2, "cm"),
                                                legend_height = unit(0.05, "mm"),
                                                values_gp     = gpar(fontsize=8),
                                                legend_direction = "horizontal") 
                  , row_names_gp = gpar(fontsize=term_font_size)
                  # , row_names_max_width = unit(5, "cm")
                  # , width = unit(2,'cm')
                  , row_names_side = "left"
                  , column_names_gp = gpar(fontsize=8)
                  , rect_gp = gpar(col = 'black')
                  , na_col = 'white'
                  # , ...
                  )
    
    ComplexHeatmap::draw(hm, heatmap_legend_side = "top")
    
  }
  
  if(is.list(top_enrichment) & !is.data.frame(top_enrichment)) {
    m.tmp <- lapply(top_enrichment, get_matrix)
    go.terms <- unique(unlist(lapply(m.tmp, rownames)))
    m <- matrix(nrow = length(go.terms)
                , ncol = length(m.tmp)
                , dimnames = list(go.terms, names(m.tmp)))
    for(i in rownames(m)) {
      for(j in colnames(m)) {
        m[i,j] <- tryCatch(-log10(m.tmp[[j]][i,]), error = function(e) return(NA))
      }
    }
  } else {
    m <- get_matrix(top_enrichment)
  }
  if(transpose) m <- t(m)
  plot.heatmap(m, ...)
}

# Modelling ====
find.variable.genes <- function(m
                                , n=1000
                                , x0 = c(-0.5, 0.5)
                                , ret.plot=T)
{
  
  log2_cv2 <- log2(apply(m, 1, function(r) (sd(r)/mean(r))**2))
  log2_mn  <- log2(apply(m, 1, function(r) mean(r)))
  
  idx      <- names(which(!is.na(log2_cv2)))
  
  log2_cv2 <- log2_cv2[idx]
  log2_mn  <- log2_mn[idx]
  
  noise.model <- function(log2_mn, log2_cv2)
  {
    function(x) sum(abs(log2((2 ** log2_mn) ** x[1] + x[2]) - log2_cv2))
  }
  
  xopt <- optim(x0, noise.model(log2_mn, log2_cv2))
  log2_cv2_fit  <- log2(( 2 ** log2_mn) ** xopt$par[1] + xopt$par[2])
  log2_cv2_diff <- log2_cv2 - log2_cv2_fit
  
  idx.ord <- order(log2_cv2_diff, decreasing = T)
  
  log2_cv2_diff[idx.ord][1:n] <- 'TRUE'
  log2_cv2_diff[idx.ord][(n+1):length(log2_cv2_diff)] <- 'FALSE'
  
  y <- cbind.data.frame(    'log2_cv2'      = log2_cv2[idx.ord]
                          , 'log2_mean'     = log2_mn[idx.ord]
                          , 'log2_cv2_fit'  = log2_cv2_fit[idx.ord]
                          , 'High'          = log2_cv2_diff[idx.ord] )
  
  genes <- rownames(y)[y$High=='TRUE']
  
  if(ret.plot) {
    
    p <- ggplot(y, aes(x=log2_mean, y=log2_cv2, col=High)) + 
      geom_point(size=1) +
      geom_line(aes(y=log2_cv2_fit), col='#CC0000') +
      theme_bw() + my_theme + scale_color_manual(values=c('black', '#0066CC')) +
      xlab("log2(mean)") + ylab("log2(cv^2)")
    
    return(list("genes" = genes,
                "plot" = p))
    
  } else {
    return(genes)
  }
}
