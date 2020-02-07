# Quality Controls - Raw data ----
read_salmon_stats <- function(DATADIR)
{
  fl <- list.files(DATADIR
                   , pattern = 'meta_info.json'
                   , recursive = T
                   , full.names = T)
  nm <- list.files(DATADIR)
  stats <- lapply(fl, function(i) {
    tmp <- rjson::fromJSON(file = i)
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
    stats$sample <- NULL
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

calc_gene_biotype_stats <- function(DATADIR, m, gene_info
                                    , type = 'salmon')
{
  
  if(missing(m)) {
    message('[+] Loading expression matrix')
    
    if(type=='salmon') {
      message(" -- salmon")
      m <- read_salmon_quant(DATADIR = DATADIR)
      
    } else if(type=='hisat') {
      message(" -- hisat")
      m <- loadCounts(countsDIR = DATADIR, pattern = ".*txt$")
    }
  }
  
  if(any(grepl('mm10', DATADIR))) {
    gene_idx <- 'ensembl_gene_id'
  } else if(any(grepl('gencode|GRCh38', DATADIR))){
    gene_idx <- 'gene_name'
    type_idx <- 'gene_type'
  } else {
    gene_idx <- 'external_gene_name'
    type_idx <- 'gene_biotype'
  }
  
  if(is.character(gene_info)) {
    message(" -- reading gene info from: ", gene_info)
    gene_info <- read.delim2(gene_info, header = T, stringsAsFactors = F)
  }
  
  # biotypes statistics ---
  # mtidx <- c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
  #            "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
  #            "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
  #            "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
  #            "ENSG00000198840")
  
  biotypes <- c('protein_coding', 'lincRNA|lnc', "rRNA")
  names(biotypes) <- c('protein_coding','lncRNA','rRNA')
  biotypes_idx    <- lapply(biotypes, function(i) intersect(gene_info[which(grepl(i,gene_info[,type_idx])), gene_idx],rownames(m)))
  biotypes_idx[['mitochondrial']] <- intersect(grep("^MT-|^mt-",gene_info$gene_name, value = T),rownames(m))
  # biotypes_idx[['mitochondrial']] <- subset(gene_info, gene_id%in%mtidx)$gene_name
  biotypes_stats  <- lapply(names(biotypes_idx), function(i) 
  {
    x <- biotypes_idx[[i]]
    df <- data.frame('perc'=round(100*colSums(m[x,])/colSums(m),2), 'nr'=colSums(m[x,]))
    colnames(df) <- paste0(colnames(df), "_",i)
    return(df)
  })
  biotypes_stats <- do.call(cbind, biotypes_stats)
  
  biotypes_stats$mapped_reads <- floor(colSums(m))
  biotypes_stats$ngenes <- colSums(apply(m, 2, function(i) i>0))
  
  return(biotypes_stats)
}

exclude_failed_stats <- function(x)
{
  which(unlist(lapply(x, function(i) length(i) == 1)))
}

plot_stats <- function(stats, vars, ptitle
                       , yintercept = 50
                       , xintercept = 1e5
                       , lim = c(0,100)
                       , pal = NULL
                       , ydirection = 'less' # greater
                       )
{
  if(missing(vars)) {
    bw <- function(x) (2 * IQR(x) / length(x)^(1/3)) # Freedman–Diaconis rule
    toplot <- reshape2::melt(stats)
    
    p <- ggplot(toplot, aes(x = value, fill = variable)) +
      geom_histogram(col='black', size = 0.25, alpha = 0.8, binwidth = bw) +
      facet_grid(~variable, scales = 'free') + theme_bw() + my_theme_2 +
      scale_fill_d3()
    
  } else {
    
    if(!missing(ptitle)) {
      stats$title <- ptitle
    }
    
    if(is.null(pal)) pal <- ggsci::pal_d3()(length(unique(stats[,vars[3]])))
    
    stats[which(is.nan(stats[,vars[1]])),vars[1]] <- 0
    stats[which(is.nan(stats[,vars[2]])),vars[2]] <- 0
    
    if(!is.null(xintercept)) {
      sel_var1 <- stats[,vars[1]]>=xintercept
    } else {
      sel_var1 <- rep(TRUE, nrow(stats))
    }
    
    if(!is.null(yintercept)) {
      if(ydirection=="greater") {
        sel_var2 <- stats[,vars[2]]>=yintercept  
      } else if(ydirection=="less") {
        sel_var2 <- stats[,vars[2]]<=yintercept
      }
    } else {
      sel_var2 <- rep(TRUE, nrow(stats))
    }
    
   
    stats$pass <- sel_var1 & sel_var2
     
    p <- ggplot(stats, aes_string(x = vars[1], y = vars[2])) +
      geom_point(aes_string(col=vars[3], alpha = "pass"), size = 1) +
      # facet_grid(~variable, scales = 'free') +
      geom_hline(yintercept = yintercept, linetype = 'dashed', size = 0.25) +
      geom_vline(xintercept = xintercept, linetype = 'dashed', size = 0.25) +
      theme_bw() + my_theme_2 + scale_alpha_manual(values = c(0.2, 0.8)) +
      scale_color_manual(values = pal) + ylim(lim) + guides(alpha = FALSE)
    
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
      nm <- gsub("path_|p_","qc_",i)
      outfile <- paste0(outdir,"/",nm,".pdf")
      message(" -- saving to file: ", outfile)
      tryCatch({pdf(file = outfile, paper = 'a4', ...)
        print(p[[i]])
        dev.off()}
        , error = function(e) message(e))
    }
  } else {
    nm <- gsub("path_|p_","qc_",names(p))
    outfile <- paste0(outdir,"/",nm,".pdf")
    pdf(file = outfile, paper = 'a4', ...)
    print(p)
    dev.off()
  }
  
}

plot_all_stats <- function(gene_stats, outdirs
                           , th_mapped_reads = 100000
                           , th_detected_genes = 2000
                           , th_protein_coding = 75
                           , th_mitochondrial = 25
                           , return_plots = T
                           , save_plots = T
                           ,...)
{
  # detected genes  ---
  p_detected <-  lapply(seq_along(gene_stats), function(i) {
    message("[+] Plotting detected genes stats from: ", names(gene_stats)[i])
    tryCatch(plot_stats(gene_stats[[i]]
                        , vars = c('mapped_reads', 'ngenes', 'condition')
                        , ptitle = 'Detected genes'
                        , yintercept = th_detected_genes
                        , xintercept = th_mapped_reads
                        , lim = c(0,max(gene_stats[[i]]$ngenes))) 
             , error = function(e) return(NA))
  })
  names(p_detected) <- names(gene_stats)
  
  # % protein coding ---
  p_protein <-  lapply(seq_along(gene_stats), function(i) {
    message("[+] Plotting protein coding gene stats from: ", names(gene_stats)[i])
    tryCatch(plot_stats(gene_stats[[i]]
                        , vars = c('ngenes', 'perc_protein_coding', 'condition')
                        , ptitle = 'Protein coding genes'
                        , yintercept = th_protein_coding
                        , xintercept = th_detected_genes
                        , lim = c(0,100)) + ylab("% of mapped reads")
             , error = function(e) return(NA))
  })
  names(p_protein) <- names(gene_stats)
  
  # % mitochondrial ---
  p_mito <-  lapply(seq_along(gene_stats), function(i) {
    message("[+] Plotting Mitochondrial gene stats from: ", names(gene_stats)[i])
    tryCatch(plot_stats(gene_stats[[i]]
                        , vars = c('ngenes', 'perc_mitochondrial', 'condition')
                        , ptitle = 'Mitochondrial genes'
                        , yintercept = th_mitochondrial
                        , xintercept = th_detected_genes
                        , lim = c(0,30)) + ylab("% of mapped reads")
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

plot_all_stats_v2 <- function(gene_stats, outdir
                           , th_mapped_reads = 100000
                           , th_detected_genes = 2000
                           , th_protein_coding = 75
                           , th_mitochondrial = 25
                           , return_plots = T
                           , save_plots = T
                           , col_by = 'condition'
                           ,...)
{
  # detected genes  ---
  p_detected <- tryCatch(plot_stats(gene_stats
                                     , vars = c('mapped_reads', 'ngenes', col_by)
                                     , ptitle = 'Detected genes'
                                     , yintercept = th_detected_genes
                                     , xintercept = th_mapped_reads
                                     , lim = c(0,max(gene_stats$ngenes))
                                     , ydirection = "greater") + 
                           ylab("n. of genes") + xlab("mapped reads") + guides(col=guide_legend(nrow=2))
                          , error = function(e) {message(e);return(NA)})
  
  # % protein coding ---
  p_protein <- tryCatch(plot_stats(gene_stats
                                    , vars = c('mapped_reads', 'perc_protein_coding', col_by)
                                    , ptitle = 'Protein coding genes'
                                    , yintercept = th_protein_coding
                                    , xintercept = th_mapped_reads
                                    , lim = c(0,100)
                                    , ydirection = "greater") + 
                          ylab("% of reads") + xlab("mapped reads") + guides(col=guide_legend(nrow=2))
                         , error = function(e) {message(e);return(NA)})
  
  # % mitochondrial ---
  p_mito <- tryCatch(plot_stats(gene_stats
                                 , vars = c('mapped_reads', 'perc_mitochondrial', col_by)
                                 , ptitle = 'Mitochondrial genes'
                                 , yintercept = th_mitochondrial
                                 , xintercept = th_mapped_reads
                                 , lim = c(0,100)
                                 , ydirection = "less") + 
                       ylab("% of reads") + xlab("mapped reads") + guides(col=guide_legend(nrow=2))
                      , error = function(e) {message(e);return(NA)})
  
  # Save QC results ---
  plots <- list('p_detected' = p_detected
                ,'p_protein' = p_protein
                ,'p_mito'    = p_mito)
  
  if(save_plots) {
    message("[+] Saving QC results ...")
    save_qc_plots(p = plots, outdir = outdir, ...)
    message(" -- done.")
  }
  
  if(return_plots) return(plots)
}

plot_hisat_stats <- function(stats, size_cutoff)
{
  toplot <- stats
  idx <- toplot[with(toplot, order(library_size, decreasing = F)),'sample']
  toplot <- reshape2::melt(toplot[,c(2:4,7)], id.var = 'sample')
  toplot$sample <- factor(toplot$sample, levels = idx)
  toplot$variable <- factor(toplot$variable, levels = c('unmapped','multi_map','uniquely_map'))
  
  p0 <- ggplot(toplot, aes(x = sample, y = value, fill = variable, col = variable)) +
    geom_col(size = 0.25, alpha = 0.9) +   
    ylab("Number of Reads") 
 
  
  if(!missing(size_cutoff)) {
    p0 <- p0 +
      geom_hline(yintercept = size_cutoff, linetype = 'dashed', col = 'red', alpha = 0.5, size=0.25)
  }
  p <- p0 +
    theme_classic() + my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank()) +
    scale_fill_d3() + scale_color_d3()
  
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
