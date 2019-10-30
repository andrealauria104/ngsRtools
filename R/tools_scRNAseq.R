# Tools for scRNA-seq data analysis
#source("tools_RNAseq.R")

require(SingleCellExperiment)
require(scater)
require(scran)
require(M3Drop)
require(SC3)
require(monocle)
require(Seurat)

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
  # miss <- setdiff(nm,sapply(strsplit(fl, "/"), "[[",13))
  # miss <- setdiff(nm,sapply(strsplit(fl, "/"), function(i) grep("EPI_|Epi_|ME_|EMPTY_|Minibulk_", i, value = T)))
  # if(length(miss)>0) {
  #   miss_idx <- which(nm%in%miss)
  #   counts[[length(fl)+1]] <- rep(0,length(counts[[1]]))
  #   names(counts) <- nm[-miss_idx]
  #   names(counts)[length(counts)] <- miss
  #
  #   counts <- counts[c((1:miss_idx-1),length(counts), (miss_idx):(length(counts)-1))]
  #
  # } else {
  #   names(counts) <- nm
  # }
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
    cds <- newCellDataSet(normcounts(sce)
                          , phenoData   = pd
                          , featureData = fd
                          , expressionFamily = VGAM::negbinomial.size())
    sizeFactors(cds) <- sizeFactors(sce)
  } else {
    # Use monocle normalization
    cds <- newCellDataSet(counts(sce)
                          , phenoData   = pd
                          , featureData = fd
                          , expressionFamily = VGAM::negbinomial.size())
    cds <- estimateSizeFactors(cds)
  }

  cds <- estimateDispersions(cds)

  return(cds)
}

# Quality Controls - Raw data ----
read_salmon_stats <- function(DATADIR)
{
  require(rjson)
  fl <- list.files(DATADIR
                   , pattern = 'meta_info.json'
                   , recursive = T
                   , full.names = T)
  nm <- list.files(DATADIR)
  stats <- lapply(fl, function(i) {
    tmp <- fromJSON(file = i)
    idx <- grep('mapped|processed', names(tmp))
    tmp <- tmp[idx]
    tmp <- cbind.data.frame(tmp)
    return(tmp)
  } )

  stats <- do.call(rbind, stats)
  tmp  <- process_barcode(nm)
  stats <- cbind.data.frame(stats, tmp)

  return(stats)
}

read_stats <- function(DATADIR, type = 'salmon')
{
  if(type=='salmon') {
    return(read_salmon_stats(DATADIR))
  } else if(type=='hisat') {
    stats   <- AlignStat(STATDIR = DATADIR)
    tmp  <- process_barcode(stats$sample)
    tmp$sample <- NULL
    stats <- cbind.data.frame(stats, tmp)
  }
}

get_gene_stats <- function(DATADIR, gene_info
                           , type = 'salmon')
{

  message('[+] Loading expression matrix')

  if(type=='salmon') {
    message(" -- salmon")
    m <- read_salmon_quant(DATADIR = DATADIR)

  } else if(type=='hisat') {
    message(" -- hisat")
    m <- loadCounts(countsDIR = DATADIR, pattern = ".*txt$")

  }

  if(any(grepl('mm10', DATADIR))) {
    gene_idx <- 'ensembl_gene_id'
  } else if(any(grepl('gencode_M1', DATADIR))){
    gene_idx <- 'gene_name'
    type_idx <- 'gene_type'
  } else {
    gene_idx <- 'external_gene_name'
    type_idx <- 'gene_biotype'
  }

  pc_idx <- gene_info[which(gene_info[,type_idx]=='protein_coding'), gene_idx]
  mt_idx <- gene_info[which(grepl('Mt_',gene_info[,type_idx]) | grepl('^mt-',gene_info[,gene_idx])), gene_idx]
  ln_idx <- gene_info[which(grepl('lincRNA|lnc',gene_info[,type_idx])), gene_idx]
  gene_stats <- list()
  gene_stats$mapped_reads <- floor(colSums(m))
  gene_stats$ngenes <- colSums(apply(m, 2, function(i) i>0))
  gene_stats$protein_coding <- round(100*colSums(apply(m[intersect(pc_idx,rownames(m)),], 2, function(i) i>0))/gene_stats$ngenes,2)
  ln <- intersect(mt_idx,rownames(m))
  if(length(ln)>0) {
    gene_stats$mitochondrial <- round(100*colSums(apply(m[ln,], 2, function(i) i>0))/gene_stats$ngenes,2)
  }
  gene_stats$noncoding <- round(100*colSums(apply(m[intersect(ln_idx,rownames(m)),], 2, function(i) i>0))/gene_stats$ngenes,2)

  gene_stats <- do.call(cbind.data.frame, gene_stats)
  tmp  <- process_barcode(rownames(gene_stats))
  gene_stats <- cbind.data.frame(gene_stats, tmp)

  return(gene_stats)
}

exclude_failed_stats <- function(x)
{
  which(unlist(lapply(x, function(i) length(i) == 1)))
}

plot_stats <- function(stats, vars, ptitle
                       , yintercept = 50
                       , xintercept = 1e5
                       , lim = c(0,100))
{
  #source("theme_setting.R")

  if(missing(vars)) {
    bw <- function(x) (2 * IQR(x) / length(x)^(1/3)) # Freedman–Diaconis rule
    toplot <- reshape2::melt(stats)

    p <- ggplot(toplot, aes(x = value, fill = variable)) +
      geom_histogram(col='black', size = 0.25, alpha = 0.8, binwidth = bw) +
      facet_grid(~variable, scales = 'free') + theme_classic() + my_theme +
      scale_fill_d3()

  } else {

    if(!missing(ptitle)) {
      stats$title <- ptitle
    }

    p <- ggplot(stats, aes_string(x = vars[1], y = vars[2])) +
      geom_point(aes_string(col=vars[3]), size = 1.5, alpha = 0.8) +
      # facet_grid(~variable, scales = 'free') +
      geom_hline(yintercept = yintercept, linetype = 'dashed', size = 0.25) +
      geom_vline(xintercept = xintercept, linetype = 'dashed', size = 0.25) +
      theme_classic() + my_theme +
      scale_color_d3() + ylim(lim)

    if(!missing(ptitle)) {
      p <- p + facet_grid(~title)
    }
  }

  return(p)
}

save_qc_plots <- function(p, outdir, ...)
{
  if(length(p)>1) {
    if(!dir.exists(outdir)) dir.create(outdir)
    for(i in names(p)) {
      i <- gsub("path_","",i)
      outfile <- paste0(outdir,"/",i,".pdf")
      tryCatch({pdf(file = outfile, paper = 'a4', ...)
        print(p[[i]])
        dev.off()}
        , error = function(e) message(e))
    }
  } else {
    nm <- gsub("path_","",names(p))
    outfile <- paste0(outdir,"_",nm,".pdf")
    pdf(file = outfile, paper = 'a4', ...)
    print(p)
    dev.off()
  }

}

plot_all_stats <- function(gene_stats, outdirs, return_plots=T, ...)
{
  # detected genes  ---
  p_detected <-  lapply(seq_along(gene_stats), function(i) {
    message("[+] Plotting detected genes stats from: ", names(gene_stats)[i])
    tryCatch(plot_stats(gene_stats[[i]]
                        , vars = c('mapped_reads', 'ngenes', 'condition')
                        , ptitle = 'Detected genes'
                        , yintercept = 2000
                        , xintercept = 100000
                        , lim = c(0,max(gene_stats[[i]]$ngenes)))
             , error = function(e) return(NA))
  })
  names(p_detected) <- names(gene_stats)

  # % protein coding ---
  p_protein <-  lapply(seq_along(gene_stats), function(i) {
    message("[+] Plotting protein coding gene stats from: ", names(gene_stats)[i])
    tryCatch(plot_stats(gene_stats[[i]]
                        , vars = c('ngenes', 'protein_coding', 'condition')
                        , ptitle = 'Protein coding genes'
                        , yintercept = 75
                        , xintercept = 2000
                        , lim = c(0,100))
             , error = function(e) return(NA))
  })
  names(p_protein) <- names(gene_stats)

  # % mitochondrial ---
  p_mito <-  lapply(seq_along(gene_stats), function(i) {
    message("[+] Plotting Mitochondrial gene stats from: ", names(gene_stats)[i])
    tryCatch(plot_stats(gene_stats[[i]]
                        , vars = c('ngenes', 'mitochondrial', 'condition')
                        , ptitle = 'Mitochondrial genes'
                        , yintercept = 10
                        , xintercept = 2000
                        , lim = c(0,30))
             , error = function(e) return(NA))
  })
  names(p_mito) <- names(gene_stats)

  # Save QC results ---
  plots <- list('p_detected' = p_detected
                ,'p_protein' = p_protein
                ,'p_mito'    = p_mito)

  message("[+] Saving QC results ...")
  mapply(save_qc_plots, p = plots, outdir = outdirs, ...)
  message(" -- done.")
  dev.off()

  if(return_plots) return(plots)
}

plot_hisat_stats <- function(stats, size_cutoff)
{
  #source("theme_setting.R")
  toplot <- stats
  idx <- toplot[with(toplot, order(library_size, decreasing = F)),'sample']
  toplot <- reshape2::melt(toplot[,c(2:4,7)], id.var = 'sample')
  toplot$sample <- factor(toplot$sample, levels = idx)
  toplot$variable <- factor(toplot$variable, levels = c('unmapped','multi_map','uniquely_map'))

  p0 <- ggplot(toplot, aes(x = sample, y = value, fill = variable)) +
    geom_col(col='black', size = 0.25, alpha = 0.8) + ylab("Number of Reads")

  if(!missing(size_cutoff)) {
    p0 <- p0 +
      geom_hline(yintercept = size_cutoff, linetype = 'dashed', col = 'red', alpha = 0.5)
  }
  p <- p0 +
    theme_light() + my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_d3()

  return(p)
}

save_hs_pstat <- function(x, outfile)
{
  pdf(file = outfile, paper = 'a4r', w=20, h=10)
  print(x)
  dev.off()
}

get_qc_cells <- function(x
                         , maprds_th = 100000
                         , ngenes_th = 2500)
{
  qc_cells <- x[with(x, mapped_reads >= maprds_th & ngenes >= ngenes_th), ]
}

# t-distributed Stochastic Neighbor Embedding ----
getTSNE <- function(x
                    , groupBy    = NULL
                    , pal        = NULL
                    , point_size = 1.5
                    , marker     = NULL
                    , ...)
{
  require(Rtsne)
  #source("theme_setting.R")

  tsne <- Rtsne(x, ...)
  tsneplot <- as.data.frame(tsne$Y)
  rownames(tsneplot) <- rownames(x)

  if( !is.null(groupBy) ) {
    if(!is.list(groupBy)) {
      tsneplot$sample <- groupBy[rownames(tsneplot)]
    } else {
      tsneplot$sample     <- groupBy[[1]][rownames(tsneplot)]
      tsneplot$additional <- groupBy[[2]][rownames(tsneplot)]
    }

  } else {
    tsneplot$sample <- rownames(tsneplot)
  }

  if(is.null(marker)) {
    p0 <- ggplot(tsneplot, aes(x=V1, y=V2, col=sample)) +
      geom_point(size=point_size) + scale_color_manual(values=pal)
  } else {
    mnames <- unique(marker$name)
    if(length(mnames)>1) {
      tmp <- list()
      for(i in mnames) {
        tmp[[i]] <- tsneplot
        tmp_mark <- subset(marker, name==i)
        idx <- match(rownames(tmp[[i]]), tmp_mark$cell)
        tmp[[i]]$marker <- i
        tmp[[i]]$expression <- tmp_mark$expression[idx]
      }
      tsneplot <- do.call(rbind, tmp)
      rm(tmp_mark)
    } else {
      idx <- match(rownames(tsneplot), marker$cell)
      tsneplot$marker <- marker$name[idx]
      tsneplot$expression <- marker$expression[idx]
    }

    p0 <- ggplot(tsneplot, aes(x=V1, y=V2, col=expression)) +
      geom_point(size=point_size) + facet_wrap(~marker, ncol = 3) +
      scale_color_gradientn(colours = RColorBrewer::brewer.pal(4,"Reds")
                            , guide = guide_colourbar(barheight = 0.6, title.vjust = 1 )
                            , values = c(0.01,0.3,0.6,1)
                            # , limits = c(0,5)
                            )
  }

  p <- p0 + xlab("t-SNE 1") + ylab("t-SNE 2") +
    theme_bw() + my_theme +
    theme(panel.grid = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1
          , strip.background = element_blank()
          , strip.text = element_text(size = 8, face = "bold"))


  return(list("tsne" = tsne, "plot" = p))

}

# Dimensionality reduction ----
plot_dimred_cells <- function(cds, pal, ...)
{
  #source("theme_setting.R")
  if(missing(pal)) pal <- 'grey'
  p <- plot_cells(cds, ...) +
    scale_color_manual(values = pal) +
    theme_bw() + my_theme + theme(panel.grid = element_blank())

  return(p)
}

# Modeling ----
model_variance <- function(sce
                           , assay.type  = "logcounts"
                           , method      = "loess"
                           , ntop        = 2500
                           , fdr.cutoff  = 0.05
                           , bio.cutoff  = 0.5
                           , mean.cutoff = 0.1)
{

  require(scran)
  plot_variance <- function(decomp)
  {
    #source("theme_setting.R")

    p <- ggplot(as.data.frame(decomp), aes(x=mean, y=total, col=hv)) +
      geom_point(size=1) +
      geom_line(aes(y=tech), col='#CC0000') +
      theme_bw() + my_theme + scale_color_manual(values=c('black', '#0066CC')) +
      xlab("Mean - log2[CPM]") + ylab("Variance")

    return(p)
  }

  # Fit a mean-dependent trend to the gene-specific variances (technical variance)
  varfit <- trendVar(sce
                     , parametric = T
                     , method     = method
                     , assay.type = assay.type
                     , use.spikes = F)

  # Decompose the gene-specific variance into biological and technical components
  decomp  <- decomposeVar(sce, varfit)
  decomp  <- decomp[order(decomp$bio, decreasing=TRUE), ]

  hv      <- decomp

  if(!is.null(mean.cutoff)) {
    hv <- subset(hv, mean >= mean.cutoff)
  }
  if(!is.null(fdr.cutoff)) {
    hv <- subset(hv, FDR <= fdr.cutoff)
  }
  if(!is.null(bio.cutoff)) {
    hv <- subset(hv, bio >= bio.cutoff)
  }

  hv      <- hv[order(hv$bio, decreasing=TRUE), ]

  if(is.null(ntop) || ntop > nrow(hv)) {
    hvgenes <- rownames(hv)
  } else {
    hvgenes <- rownames(hv)[1:ntop]
  }

  decomp$hv <- F
  decomp$hv[match(hvgenes, rownames(decomp))] <- T

  p <- plot_variance(decomp)

  return(list("model" = decomp, "hvgenes" = hvgenes, "plot" = p))
}

regress_variation_out <- function(expression_matrix, covariates, model_formula)
{
  mod <- model.matrix(as.formula(model_formula)
                      , data = covariates
                      , drop.unused.levels = TRUE)
  fit <- limma::lmFit(expression_matrix, mod)
  beta <- fit$coefficients[, -1, drop = FALSE]
  beta[is.na(beta)] <- 0
  expression_matrix <- as.matrix(expression_matrix) - beta %*% t(mod[, -1])

  return(expression_matrix)
}
# Differential Expression ----
plot_cmarker_expression <- function(ebs, genes, pal, structure = "Cluster",log = T, assay.type = NULL, point.size = 0.5, plot.type = 'boxplot', filter = F)
{
  #source("theme_setting.R")

  if(is.null(assay.type)) {
    if(log) {
      marker <- logcounts(ebs)[genes,]
    } else {
      marker <- normcounts(ebs)[genes,]
    }
  } else {
    marker <- assay(ebs, assay.type)[genes,]
  }

  marker <- reshape2::melt(marker)
  colnames(marker) <- c("name","cell","expression")
  marker$cluster <- colData(ebs)[marker$cell,grep(structure, colnames(colData(ebs)))]
  if(filter) marker <- marker[marker$expression>0,]
  if(plot.type=='boxplot') {
    dodge <- position_dodge(width = 0.8)
    p0 <- ggplot(marker, aes(x=cluster, y=expression, col = cluster)) + geom_boxplot(lwd = 0.25, outlier.size = 0.1) + geom_jitter(size = point.size)
  } else if(plot.type=='violin') {
    dodge <- position_dodge(width = 0.8)
    p0 <- ggplot(marker, aes(x=cluster, y=expression, col = cluster)) +
      geom_violin(trim = T, position = dodge) + geom_boxplot(width=0.08, outlier.color = NA, position = dodge) +
      stat_summary(fun.y=median, geom="point", size=0.5, color="black", position = dodge)
  }
  p <- p0 +
    facet_wrap(~name, ncol = 3) +
    theme_bw() + my_theme + scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) + xlab(structure) + ylab("Expression")  +
    theme(panel.grid = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1
          , strip.background = element_blank()
          , strip.text = element_text(size = 8, face = "bold"))

  return(p)
}

run_differrential_expression <- function(i, cds, structure, ntop, direction = "up", rank.var = "pval", ...)
{
  pData(cds)$category <- pData(cds)[,structure] == i
  cgroup_i <- rownames(pData(cds))[which(pData(cds)$category==T)]
  cgroup_rest <- rownames(pData(cds))[which(pData(cds)$category==F)]
  mean_1 <- rowMeans(exprs(cds)[,cgroup_i])
  mean_2 <- rowMeans(exprs(cds)[,cgroup_rest])

  test  <- differentialGeneTest(cds, fullModelFormulaStr = "~ category", ...)
  test$mean_cgroup_i <- mean_1[rownames(test)]
  test$mean_cgroup_rest <- mean_2[rownames(test)]
  test$logFC <- log2(test$mean_cgroup_i+1) - log2(test$mean_cgroup_rest+1)

  if(direction=="up") {
    results <- subset(test, mean_cgroup_i>mean_cgroup_rest)
  } else if(direction=="dw") {
    results <- subset(test, mean_cgroup_i<mean_cgroup_rest)
  } else if(direction=="all") {
    results <- test
  }
  if(!is.null(ntop)) results <- head(results[order(results[,rank.var], decreasing = rank.var=="logFC"),], n = ntop)

  return(results)
}

find_marker_genes <- function(cds, structure = "Cluster", ntop = 50, names_only = T, ...)
{
  structures <- levels(pData(cds)[,structure])
  structure_markers <- lapply(structures, run_differrential_expression
                              , cds  = cds
                              , ntop = ntop
                              , structure = structure
                              , ...)
  names(structure_markers) <- structures

  if((!is.null(ntop) || ntop > 0) && names_only == T) {
    top <- lapply(structure_markers, "[[", "feature_symbol")
    top <- as.data.frame(do.call(cbind, top))

    return(top)
  } else {
    return(structure_markers)
  }


}
# Pseudotime Analysis ----
plot_pseudotemporal_ordering <- function(cds, ffeature, cfeature, pal, ...)
{
  #source("theme_setting.R")

  if(missing(pal)) pal <- ggsci::pal_aaas()(10)

  if(missing(ffeature)) {
    ffeature <- ""
  } else if(missing(cfeature)) {
    cfeature <- ffeature
  } else if(missing(cfeature) & missing(ffeature)) {
    stop(message("[!] Please, provide color/facet feature."))
  }

  if(length(ffeature)!=1 || length(cfeature)!=1) {

    if(length(cfeature)>1) {
      mergefeature <- cfeature
      cfeature <- "sample"
    } else if(length(ffeature)>1) {
      mergefeature <- ffeature
      ffeature <- "sample"
    }

    nm <- apply(phenoData(cds)@data[,mergefeature], 1, paste, collapse = "_")
    phenoData(cds)$sample <- gsub("_","/",nm)

    p0 <- plot_cell_trajectory(cds, color_by = cfeature, ...) +
      facet_wrap(paste0("~",ffeature), nrow = 1)

  } else {
    p0 <- plot_cell_trajectory(cds, color_by = cfeature, ...)
  }

  if(is.numeric(phenoData(cds)@data[,cfeature])) {

    p0 <- p0 + scale_color_gradientn(colours = RColorBrewer::brewer.pal(4,"Purples")
                                     , guide = guide_colourbar(barheight = 0.6, title.vjust = 1 )
                                     # , values = c(0.01,0.3,0.6,1)
                                     # , limits = c(-1,1)
    )
  } else {
    p0 <- p0 + scale_color_manual(values = pal)
  }
  pct <- p0 + theme_bw() + my_theme +
    theme(panel.grid = element_blank()
          , strip.background = element_blank()
          , strip.text = element_text(size = 8, face = "bold"))

  return(pct)
}

plot_cells_in_pseudotime <- function(cds, structure, pal)
{
  #source("theme_setting.R")

  if(length(structure)>1) {
    nm <- apply(pData(cds)[,structure], 1, paste, collapse = "_")
    pData(cds)$sample <- gsub("_","/",nm)
    structure   <- "sample"
  }

  ptime  <- pData(cds)[,c("Pseudotime",structure)]
  pptime <- ggplot(ptime, aes_string(x="Pseudotime", y=structure, col=structure)) + geom_point() +
    theme_bw() + my_theme + scale_color_manual(values = pal)

  return(pptime)
}

analyze_beam_cluster_profile <- function(pbeam)
{
  #source("theme_setting.R")
  x <- pbeam$heatmap_matrix
  colnames(x) <- c(-99:100)
  beam_clustering <- pbeam$annotation_row
  cp <- lapply(unique(beam_clustering$Cluster),
               function(i)
               {
                 y <- x[rownames(beam_clustering)[beam_clustering$Cluster==i],]
                 z <- apply(y, 2, summary)
                 z <- as.data.frame(t(z))
                 z$time <- as.numeric(rownames(z))
                 colnames(z)[c(2,5)] <- c("qu1","qu3")
                 z$cluster <- i
                 return(list("z"=z,"y"=y))
               } )
  cg <- lapply(cp, function(i) rownames(i$y))
  names(cg) <- paste0("cluster_",unique(beam_clustering$Cluster))
  cp <- lapply(cp, "[[", "z")
  cp <- do.call(rbind.data.frame, cp)

  p <- ggplot(cp, aes(x=time, y=Mean, col=cluster)) + geom_line() +
    facet_wrap(~cluster, ncol = 1) + xlab("Pseudotime") + geom_vline(xintercept = 0, lwd=0.25, linetype="dashed") +
    ylab("Smoothed expression") + scale_x_continuous(breaks = 0) +
    geom_ribbon(aes(ymin=qu1, ymax=qu3, fill=cluster), alpha = 0.2) +
    theme_bw() + my_theme + theme(strip.background = element_blank(), strip.text = element_blank())

  return(list("p"=p,"gene_clusters"=cg))
}
# Clustering ----
plot_cinfo <- function(cinfo, structure, feature, pal, info_type = "absolute", mode = "bar")
{
  #source("theme_setting.R")
  if(length(feature)>1) {
    nm <- apply(cinfo[,feature], 1, paste, collapse = "_")
    cinfo$sample <- gsub("_","/",nm)
    feature   <- "sample"
  }

  if(info_type == "absolute") {
    p <- ggplot(cinfo, aes_string(x=structure, fill=feature)) + geom_bar() +
      theme_bw() + my_theme +
      # theme(legend.title = element_blank()) +
      scale_fill_manual(values = pal) + xlab(structure) + ylab("n. of cells")
  } else if(info_type == "relative") {
    print("here")
    cinfo     <- split(cinfo, cinfo[,feature])
    cinfo     <- lapply(cinfo, function(x)
    {
      x$tot_n_cells = nrow(x)
      return(x)
    } )
    cinfo <- do.call(rbind, cinfo)

    cinfo     <- split(cinfo, cinfo[,c(structure,feature)])
    cinfo     <- lapply(cinfo, function(x)
    {
      x$relative_n_cells = nrow(x)/x$tot_n_cells
      return(x)
    } )

    cinfo     <- do.call(rbind, cinfo)

    if(mode=="reverse_bar") {
      p0 <- ggplot(cinfo, aes_string(x=feature, y = "relative_n_cells", fill=structure)) +
        geom_bar(position = "fill",stat = "identity")
      x_lab <- feature
    } else if(mode=="pie") {
      require(grid)
      mp <- cinfo
      mp$value <- mp$relative_n_cells
      mp <- unique(mp[,c(structure,"sample","value")])
      mp$ypos <- cumsum(mp$value)[4] - cumsum(mp$value) + mp$value/2
      mp$perc <- round(mp$value/sum(mp$value),4)
      pie <- vector(mode = 'list')

      for(s in unique(mp[,structure])) {

        mps <- mp[mp[,structure]==s,]
        mps$ypos <- cumsum(mps$value)[nrow(mps)] - cumsum(mps$value) + mps$value/2
        mps$perc <- round(mps$value/sum(mps$value),4)

        pie[[s]] <- ggplot(mps, aes(x="", y=value, fill=sample)) +
          geom_bar(width = 1, stat = "identity", col = 'black') +
          coord_polar("y", start=0)  + scale_fill_manual(values = pal) +
          blank_theme + ggtitle(s) +
          theme(axis.text.x=element_blank()
                , plot.title = element_text(size = 8, face = "plain", hjust = 0.5, vjust = 0)
                , text = element_text(size = 8)
                , legend.position = 'none'
                , plot.margin=unit(c(0,0,0,0),"cm")
                , panel.spacing = unit(c(0, 0, 0, 0), "cm"),)
      }
      # dev.off()
      if(structure=="State") {
        pushViewport(viewport(layout = grid.layout(4, 3, widths = unit(2,"cm"), heights = unit(2,"cm"))))
        vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
        grid.text("Relative proportion of cells by state", vp = vplayout(1, 2), gp = gpar(fontsize = 8))
        print(pie$A, vp = vplayout(2, 1))
        print(pie$E, vp = vplayout(2, 3))
        print(pie$C, vp = vplayout(3, 2))
        print(pie$B, vp = vplayout(4, 1))
        print(pie$D, vp = vplayout(4, 3))
      } else if(structure=="Cluster") {
        pushViewport(viewport(layout = grid.layout(4, 3, widths = unit(2,"cm"), heights = unit(2,"cm"))))
        vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
        grid.text("Relative proportion of cells by cluster", vp = vplayout(1, 2), gp = gpar(fontsize = 8))
        print(pie[[1]], vp = vplayout(2, 1))
        print(pie[[2]], vp = vplayout(2, 2))
        print(pie[[3]], vp = vplayout(2, 3))
        print(pie[[6]], vp = vplayout(3, 1))
        print(pie[[4]], vp = vplayout(3, 2))
        print(pie[[5]], vp = vplayout(3, 3))
      }

      # gridExtra::grid.arrange(pie$A,pie$B,pie$C,pie$D,pie$E, nrow = 2)
      p0 <- NULL
    } else if(mode=="bar") {
      mp   <- cinfo
      mp   <- unique(mp[,c(structure,"sample","relative_n_cells")])
      refs <- unique(cinfo$sample)

      mp <- split(mp, mp[,structure])
      mp <- lapply(mp, function(x) {
        if(length(setdiff(refs, x$sample))>0) {
          y <- data.frame("sample" = setdiff(refs, x$sample)
                          , structure = x[,structure]
                          , "relative_n_cells" = 0)
          colnames(y)[colnames(y)=="structure"] <- structure
          x <- rbind.data.frame(x, y)
          rownames(x) <- NULL
        }
        return(x)
      })

      mp <- do.call(rbind, mp)
      mp[,structure] <- factor(mp[,structure], levels = levels(mp[,structure])[order(levels(mp[,structure]), decreasing = F)])
      p0 <- ggplot(mp, aes_string(x=structure, y = "relative_n_cells", fill=feature)) +
        geom_col(position = 'dodge')
      # + coord_flip()
      x_lab <- structure
    }

    if(mode!="pie") {
      p <- p0 +
        theme_bw() + my_theme +
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        scale_fill_manual(values = pal) + xlab(x_lab) + ylab("Percent of cells")
    } else {
      p <- p0
    }
  } else {
    stop(message("[!] Invalid info_type (absolute/relative)."))
  }

  return(p)
}

get_cell_clusters_heatmap <- function()
{

}
