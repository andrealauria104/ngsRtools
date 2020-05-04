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

process_barcode <- function(cell_barcode
                            , barcode_info # named numeric vector indexing metadata in barcode
                            , additional_info # list/vector with additional information
                            , split_char = "\\_"
                            )
{
  metadata <- list()
  metadata[['Sample']] <- cell_barcode
  
  for(i in names(barcode_info)) {
    metadata[[i]] <- sapply(strsplit(cell_barcode, split = split_char), "[[", barcode_info[i])
  }
  
  for(i in names(additional_info)) {
    if(length(additional_info[[i]]==1)) {
      metadata[[i]] <- rep(additional_info[[i]], length(cell_barcode))
    } else {
      metadata[[i]] <- additional_info[[i]]
    }
  }
  metadata <- do.call(cbind.data.frame, metadata)
  metadata$Library <- 'single_cell'
  
  # handle EMPTY/Minibulk
  if(any(grepl("empty", metadata$Sample, ignore.case = T, perl = T))) {
    empty_idx <- grep("empty", metadata$Sample, ignore.case = T, perl = T)
    metadata$Library[empty_idx] <- 'Empty_well'
  }
  if(any(grepl("(mini)?bulk", metadata$Sample, ignore.case = T, perl = T))) {
    bulk_idx <- grep("(mini)?bulk", metadata$Sample, ignore.case = T, perl = T)
    metadata$Library[bulk_idx] <- 'Minibulk'
  }
  return(metadata)
}

read_metadata <- function(metadata
                          , metadata_from_barcode = F
                          , cell_barcode
                          , barcode_info
                          , additional_info)
{
  if(is.character(metadata) && !metadata_from_barcode) {
    mdata <- read.delim(metadata, header = T, stringsAsFactors = F)
  } else if(!is.data.frame(metadata) && metadata_from_barcode) {
    mdata <- process_barcode(cell_barcode, barcode_info, additional_info)
  } else if(is.data.frame(metadata)) {
    mdata <- metadata
  }else {
    stop(message('[!] Please, provide valid path to sample metadata or process barcode.'))
  }
  return(mdata)
}

read_stats <- function(statdir
                       , metadata = NULL
                       , type = 'hisat'
                       , metadata_from_barcode = F
                       , barcode_info
                       , additional_info)
{
  if(any(grepl("multiqc",statdir)) & type == 'hisat'){
    
    if(file.exists(statdir) && !dir.exists(statdir)) {
      mstats <- read.delim(statdir, stringsAsFactors = F)
    } else if(dir.exists(statdir)) {
      stats_file <- list.files(statdir, pattern = 'bowtie2', full.names = T)
      mstats <- read.delim(stats_file, stringsAsFactors = F)
    } else {
      stop(message('[!] Please, provide valid path to alignment statistics.'))
    }
    
    mstats$Sample <- gsub(".hs2","",mstats$Sample)
    mstats$uniq_mapping_rate <- round(100*mstats$unpaired_aligned_one/mstats$total_reads,2)
    colnames(mstats) <- c('Sample','Unmapped','Uniquely_mapped'
                          ,'Multimapped','Unpaired_total','Mapping_rate'
                          ,'Total_reads', 'Uniquely_mapped_rate')
    
    mdata <- read_metadata(metadata = metadata
                           , metadata_from_barcode = metadata_from_barcode
                           , cell_barcode = mstats$Sample
                           , barcode_info = barcode_info
                           , additional_info = additional_info)
    
    stats <- cbind.data.frame(mdata, mstats[match(mdata$Sample, mstats$Sample),-1])
  } else if(type=='salmon') {
    stats <- read_salmon_stats(statdir)
  } else if(type=='hisat') {
    message('[!] Old stats summary for compatibility only, to be deprecated')
    # for compatibility, to be deprecated
    stats <- AlignStat(STATDIR = statdir)
    tmp  <- process_barcode(stats$sample)
    stats$sample <- NULL
    stats <- cbind.data.frame(stats, tmp)
  } else {
    stop(message('[!] Invalid input.'))
  }
  return(stats)
}

calc_gene_biotype_stats <- function(count_matrix, DATADIR, gene_info
                                    , type = 'hisat'
                                    , gene_idx = 'gene_name'
                                    , type_idx = 'gene_type'
                                    , metadata = NULL
                                    , metadata_from_barcode = F
                                    , barcode_info
                                    , additional_info)
{
  
  if(missing(count_matrix)) {
    message('[+] Loading expression matrix')
    
    if(type=='salmon') {
      message(" -- salmon")
      count_matrix <- read_salmon_quant(DATADIR = DATADIR)
      
    } else if(type=='hisat') {
      message(" -- hisat")
      count_matrix <- loadCounts(countsDIR = DATADIR, pattern = ".*txt$")
    }
  } else if(is.character(count_matrix) && file.exists(count_matrix)) {
    sep <- ifelse(grepl(".csv",count_matrix), ",","\t")
    count_matrix <- read.delim(count_matrix, sep = sep, stringsAsFactors = F, header = T, row.names = 1)
    count_matrix <- as.matrix(count_matrix[,-1])
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
  biotypes_idx    <- lapply(biotypes, function(i) intersect(gene_info[which(grepl(i,gene_info[,type_idx])), gene_idx],rownames(count_matrix)))
  biotypes_idx[['mitochondrial']] <- intersect(grep("^MT-|^mt-",gene_info$gene_name, value = T),rownames(count_matrix))
  # biotypes_idx[['mitochondrial']] <- subset(gene_info, gene_id%in%mtidx)$gene_name
  biotypes_stats  <- lapply(names(biotypes_idx), function(i) 
  {
    x <- biotypes_idx[[i]]
    df <- data.frame('perc'=round(100*colSums(count_matrix[x,])/colSums(count_matrix),2), 'nr'=colSums(count_matrix[x,]))
    colnames(df) <- paste0(colnames(df), "_",i)
    return(df)
  })
  biotypes_stats <- do.call(cbind, biotypes_stats)
  
  biotypes_stats$assigned_reads <- floor(colSums(count_matrix))
  biotypes_stats$ngenes <- colSums(apply(count_matrix, 2, function(i) i>0))
  
  mdata <- read_metadata(metadata = metadata
                         , metadata_from_barcode = metadata_from_barcode
                         , cell_barcode = rownames(biotypes_stats)
                         , barcode_info = barcode_info
                         , additional_info = additional_info)
  
  biotypes_stats <- cbind.data.frame(mdata, biotypes_stats[match(mdata$Sample, rownames(biotypes_stats)),])
  
  return(biotypes_stats)
}

plot_stats <- function(stats, vars, ptitle # vars = c(x-axis, y-axis, color)
                       , yintercept = 50
                       , xintercept = 1e5
                       , lim = c(0,100)
                       , pal = NULL
                       , ydirection = 'less' # greater
                       , plot.type = 'scatter')
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
    
    if(any(is.nan(stats[,vars[1]]))) stats[which(is.nan(stats[,vars[1]])),vars[1]] <- 0
    if(any(is.nan(stats[,vars[2]]))) stats[which(is.nan(stats[,vars[2]])),vars[2]] <- 0
    
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
    
    if(length(unique(stats$pass))>1) {
      alpha_values <-  c(.2, .8)
    } else {
      alpha_values <- .8
    }
    
    if(vars[3]=='') {
      stats$experiment <- "scRNA-seq"
      vars[3] <- 'experiment'
    }
    
    if(is.null(pal)) pal <- ggsci::pal_d3()(length(unique(stats[,vars[3]])))
    
    if(plot.type=="scatter") {
      p <- ggplot(stats, aes_string(x = vars[1], y = vars[2])) +
        geom_point(aes_string(col=vars[3], alpha = "pass"), size = 1) +
        # facet_grid(~variable, scales = 'free') +
        xlab(gsub("rate","rate [%]",gsub("\\_"," ",vars[1]))) + 
        ylab(gsub("rate","rate [%]",gsub("\\_"," ",vars[2]))) +
        theme_bw() + my_theme_2 + scale_alpha_manual(values = alpha_values) +
        scale_color_manual(values = pal) + ylim(lim) + 
        guides(alpha = FALSE, col=guide_legend(title = gsub("\\_"," ",vars[3])))
      if(!is.null(xintercept)) p <- p + geom_vline(xintercept = xintercept, linetype = 'dashed', size = 0.25)
      if(!is.null(yintercept)) p <- p + geom_hline(yintercept = yintercept, linetype = 'dashed', size = 0.25)
    } else if(plot.type=="boxplot") {
      p <- ggplot(stats, aes_string(x = vars[1], y = vars[2], col=vars[3])) +
        geom_boxplot(lwd=0.25, alpha=.8) + geom_jitter(size=1, width = 0.3) + 
        xlab(gsub("rate","rate [%]",gsub("\\_"," ",vars[1]))) + 
        ylab(gsub("rate","rate [%]",gsub("\\_"," ",vars[2]))) +
        theme_bw() + my_theme_2 + scale_alpha_manual(values = alpha_values) +
        scale_color_manual(values = pal) + ylim(lim) + 
        guides(alpha = FALSE, col=guide_legend(title = gsub("\\_"," ",vars[3])))
    }
    
    if(!missing(ptitle)) {
      p <- p + facet_grid(~title)
    }
  }
  
  return(p)
}

plot_mapping_stats <- function(mapping_stats, color_by
                               , xintercept = NULL
                               , yintercept = NULL
                               , p_boxplot = T
                               , ...) 
{
  p_mapping_stats <- plot_stats(mapping_stats
                                , vars = c('Total_reads', 'Mapping_rate', color_by) # x-axis, y-axis, color
                                , xintercept = xintercept # qc pass cutoff x
                                , yintercept = yintercept # qc pass cutoff y
  , ...) 
  
  p_mapping_stats_uniq <- plot_stats(mapping_stats
                                     , vars = c('Total_reads', 'Uniquely_mapped_rate', color_by)  # x-axis, y-axis, color
                                     , xintercept = xintercept # qc pass cutoff x
                                     , yintercept = yintercept # qc pass cutoff y
  , ...)
  
  if(p_boxplot) {
    p_mapping_stats_2 <- plot_stats(mapping_stats
                                    , vars = c(color_by, 'Mapping_rate', color_by) # x-axis, y-axis, color
                                    , xintercept = xintercept # qc pass cutoff x
                                    , yintercept = yintercept # qc pass cutoff y
                                    , plot.type = 'boxplot'
                                    , ...) 
    
    p_mapping_stats_uniq_2 <- plot_stats(mapping_stats
                                         , vars = c(color_by, 'Uniquely_mapped_rate', color_by)  # x-axis, y-axis, color
                                         , xintercept = xintercept # qc pass cutoff x
                                         , yintercept = yintercept # qc pass cutoff y
                                         , plot.type = 'boxplot'
                                         , ...)
    
    p_mapping_stats_arranged <- ggpubr::ggarrange(p_mapping_stats, p_mapping_stats_uniq
                                                  , p_mapping_stats_2, p_mapping_stats_uniq_2
                                                  , ncol=2, nrow=2
                                                  , common.legend = TRUE
                                                  , legend="bottom")
  } else {
    p_mapping_stats_arranged <- ggpubr::ggarrange(p_mapping_stats, p_mapping_stats_uniq
                                                  , ncol=2, nrow=1
                                                  , common.legend = TRUE
                                                  , legend="bottom")
  }
  
  
  
  return(p_mapping_stats_arranged)
}

plot_overall_alignment_stats <- function(mapping_stats, size_cutoff)
{
  custom_theme <- theme_light() + theme(text = element_text(size = 8)
                                        , line = element_line(size=0.25)
                                        , legend.key.size = unit(0.3,'cm')
                                        , legend.spacing.x = unit(0.2, 'cm')
                                        , legend.title = element_blank()
                                        , plot.title = element_text(hjust = 0.5, size = 8)
                                        , legend.position = "bottom")
  
  plot_vars <- c('Unmapped','Multimapped','Uniquely_mapped')
  idx <- mapping_stats[with(mapping_stats, order(Total_reads, decreasing = F)),'Sample']
  cidx <- grep(paste0(c("Sample",plot_vars),collapse = "|"),colnames(mapping_stats), value = T)
  cidx <- cidx[grep("rate",cidx, invert = T)]
  
  toplot <- reshape2::melt(mapping_stats[,cidx], id.var = 'Sample')
  toplot$Sample <- factor(toplot$Sample, levels = idx)
  toplot$variable <- factor(toplot$variable, levels = plot_vars)
  levels(toplot$variable) <- gsub("\\_"," ",levels(toplot$variable))
  p <- ggplot(toplot, aes(x = Sample, y = value, fill = variable, col = variable)) +
    geom_col(size = 0.25, alpha = 0.9) + ylab("Number of reads") +
    scale_fill_d3() + scale_color_d3() +
    ggtitle("Total reads") + custom_theme + theme(axis.text.x = element_blank()
                                                  , axis.ticks.x = element_blank()
                                                  ,panel.grid.major.x = element_blank()) +
    guides(col = guide_legend(reverse = TRUE, ), fill = guide_legend(reverse = TRUE))
  
  
  if(!missing(size_cutoff)) {
    p <- p +
      geom_hline(yintercept = size_cutoff, linetype = 'dashed', col = 'red', alpha = 0.5, size=0.25)
  }
  
  return(p)
}

plot_all_stats_v2 <- function(gene_stats, outdir
                              , th_assigned_reads = 100000
                              , th_detected_genes = 2000
                              , th_protein_coding = 75
                              , th_mitochondrial = 25
                              , return_plots = T
                              , save_plots = T
                              , col_by = ''
                              , pal = NULL
                              , guide_legend_nrow = 2
                              , ...)
{
  if(col_by!='') {
    p_guide_legend <- guides(col=guide_legend(nrow=guide_legend_nrow))
  } else {
    p_guide_legend <- guides(col=FALSE)
  }
  
  # detected genes  ---
  p_detected <- tryCatch(plot_stats(gene_stats
                                    , vars = c('assigned_reads', 'ngenes', col_by)
                                    , ptitle = 'Detected genes'
                                    , yintercept = th_detected_genes
                                    , xintercept = th_assigned_reads
                                    , lim = c(0,max(gene_stats$ngenes))
                                    , ydirection = "greater"
                                    , pal = pal) + 
                           ylab("n. of genes") + xlab("Assigned reads") + p_guide_legend
                         , error = function(e) {message(e);return(NA)})
  
  # % protein coding ---
  p_protein <- tryCatch(plot_stats(gene_stats
                                   , vars = c('assigned_reads', 'perc_protein_coding', col_by)
                                   , ptitle = 'Protein coding genes'
                                   , yintercept = th_protein_coding
                                   , xintercept = th_assigned_reads
                                   , lim = c(0,100)
                                   , ydirection = "greater"
                                   , pal = pal) + 
                          ylab("% of reads") + xlab("Assigned reads") + p_guide_legend
                        , error = function(e) {message(e);return(NA)})
  
  # % mitochondrial ---
  p_mito <- tryCatch(plot_stats(gene_stats
                                , vars = c('assigned_reads', 'perc_mitochondrial', col_by)
                                , ptitle = 'Mitochondrial genes'
                                , yintercept = th_mitochondrial
                                , xintercept = th_assigned_reads
                                , lim = c(0,100)
                                , ydirection = "less"
                                , pal = pal) + 
                       ylab("% of reads") + xlab("Assigned reads") + p_guide_legend
                     , error = function(e) {message(e);return(NA)})
  
  # Save QC results ---
  plots <- list('p_detected' = p_detected
                ,'p_protein' = p_protein
                ,'p_mito'    = p_mito)
  if(col_by!='') names(plots) <- paste0(names(plots),"_",col_by)
  if(save_plots) {
    message("[+] Saving QC results ...")
    save_qc_plots(p = plots, outdir = outdir, ...)
    message(" -- done.")
  }
  
  if(return_plots) return(plots)
}

# To be deprecated ----
# to be deprecated 
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
  gene_stats$assigned_reads <- floor(colSums(m))
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

# to be deprecated
exclude_failed_stats <- function(x)
{
  which(unlist(lapply(x, function(i) length(i) == 1)))
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

# to be deprecated
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

# to be deprecated
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

# to be deprecated
save_hs_pstat <- function(x, outfile)
{
  pdf(file = outfile, paper = 'a4r', w=20, h=10)
  print(x)
  dev.off()
}

# to be deprecated
get_qc_cells <- function(x
                         , maprds_th = 100000
                         , ngenes_th = 2500)
{
  qc_cells <- x[with(x, mapped_reads >= maprds_th & ngenes >= ngenes_th), ]
}
