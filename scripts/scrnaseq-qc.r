#!/usr/bin/env Rscript

# # # # # # # # # # # # # # # # # # # # # # # # # 
# single cell RNA-sequencing - quality controls #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Command line/Env variables ---
suppressWarnings(suppressMessages(require(docopt)))
'Usage:
   scrnaseq-qc.r [-q <qcmatrix> -m <metadata> -b <biotype> -c <colby> -o <outdir>]

Options:
   -q, --qcmatrix Path to QC matrix. 
   -m, --metadata Path to experiment metadata.
   -b, --biotype Gene biotype (comma separated) for diagnostic plots. 
                 Protein coding, lncRNA, rRNA and mitochondrial genes 
                 are reported by default.
   -c, --colby Color cells by feature in metadata [default: cell_type]
   -o, --outdir Output directory. [default: .]
' -> doc

opts <- docopt(doc)
if(is.null(opts$qcmatrix) || is.null(opts$metadata)) {
  message("\n[!] Missing required arguments.\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

qcmatrix <- as.character(opts$qcmatrix)
metadata <- as.character(opts$metadata)
biotype <- opts$biotype
colby   <- as.character(opts$colby)
outdir  <- as.character(opts$outdir)

# 0. Resources ----
suppressWarnings(suppressMessages(require(ggplot2)))
suppressWarnings(suppressMessages(require(ggsci)))
suppressWarnings(suppressMessages(require(plyr)))


custom_theme <- theme_light() + theme(text = element_text(size = 8)
        , line = element_line(size=0.25)
        , legend.key.size = unit(0.3,'cm')
	, legend.spacing.x = unit(0.2, 'cm')
        , legend.title = element_blank()
        , plot.title = element_text(hjust = 0.5, size = 8)
        , legend.position = "bottom")

# Read data ---
read_stats <- function(qcmatrix, metadata) 
{
  qcmat <- read.delim(qcmatrix, header = T)
  mdata <- read.delim(metadata, header = T)
  
  stats <- cbind.data.frame(mdata, qcmat[match(mdata$Sample, qcmat$Sample),-1])
  
  return(stats)
}
# Alignment statistics ---
plot_alignment_stats <- function(stats) 
{
  
  idx <- stats[with(stats, order(tot, decreasing = F)),'Sample']
  
  qc_metrics <- c('non_ribo_uniq_map','Assigned','tot_unassigned','tot_multimap', 'ribo','non_ribo_multi_map','unmap')
  names(qc_metrics) <- c(1,2,2,1,2,2,1)
  qc_metrics_fix <- c('Uniquely mapped','Assigned','Unassigned','Multimapped','Ribosomal','Non-ribo-multimap','Unmapped')
  names(qc_metrics_fix) <- c(1,2,2,1,2,2,1)
  
  pals <- list(pal_d3(), pal_jama())
  names(pals) <- 1:2
  titles <- c("Total reads","Mapped reads")
  names(titles) <- 1:2
  
  stats$tot_unassigned <- stats$Unassigned_Ambiguity+stats$Unassigned_NoFeatures
  stats$tot_multimap <- stats$ribo+stats$non_ribo_multi_map
  
  pstats <- lapply(unique(names(qc_metrics)), function(x) 
  {
    metrics <- qc_metrics[names(qc_metrics)==x]
    metrics_fix <- qc_metrics_fix[names(qc_metrics_fix)==x]
    
    toplot <- reshape2::melt(stats[,c('Sample', metrics)], id.var = 'Sample')
    tot_idx <- ddply(toplot, .(Sample), summarize, tot = sum(value))
    idx <- tot_idx[with(tot_idx, order(tot, decreasing = F)),'Sample']
    
    toplot[,'Sample'] <- factor(toplot[,'Sample'], levels = as.character(idx))
    toplot$variable <- factor(toplot$variable, levels = rev(metrics))
    levels(toplot$variable) <- rev(metrics_fix)
    
    pal <- pals[[x]](length(metrics))
    
    p <- ggplot(toplot, aes(x = Sample, y = value, fill = variable, col = variable)) +
      geom_col(size = 0.25, alpha = 0.9) + ylab("Number of reads") +
      scale_fill_manual(values = pal) + scale_color_manual(values = pal) + 
      ggtitle(titles[x]) + custom_theme + theme(axis.text.x = element_blank()
                                                , axis.ticks.x = element_blank()
                                                ,panel.grid.major.x = element_blank()) +
      guides(col = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE))
  })
  
  return(pstats)
}
# Biotype statistics ---
plot_biotype_stats <- function(stats, colby, biotype = NULL, pal = NULL)
{
  if(is.null(pal)) {
    pal <- ggsci::pal_d3()(length(unique(stats[,colby])))
    names(pal) <- unique(stats[,colby])
  }
  if(is.null(biotype)) {
    biotype <- c('protein_coding','lncRNA','mitochondrial','rRNA') 
  } else if(grepl("\\,",biotype)) {
    biotype <- unlist(strsplit(biotype,"\\,"))
  }
  # tot mitochondrial
  stats$mitochondrial <- rowSums(stats[,grep('Mt_',colnames(stats))])
  
  # % reads per biotype
  perc_biotype <- lapply(biotype, function(x) 
    {
      y <- round(100*stats[,x]/stats[,'Assigned'],2)
      y[which(is.nan(y))] <- 0
      return(y)} 
    )
  names(perc_biotype) <- paste0("perc_",biotype)
  
  pbiotypes <- lapply(names(perc_biotype), function(b) 
    {
    toplot <- cbind.data.frame(stats[,c('Sample','cell_type', 'time','Assigned')], perc_biotype[[b]])
    colnames(toplot)[5] <- b
    
    p <- ggplot(toplot, aes_string(x = 'Assigned', y = b)) +
      geom_point(aes_string(col=colby), size = 1, alpha = 0.8) + ylab("% of reads") + 
      xlab("Assigned reads") + custom_theme + scale_color_manual(values = pal) +
      ggtitle(gsub("_|perc"," ",b))  + theme(panel.grid.minor.x = element_blank())
    
  }) 
  names(pbiotypes) <- names(perc_biotype)
  
  return(pbiotypes)
}

# 1. Read data ----
stats <- read_stats(qcmatrix, metadata)

if(!dir.exists(outdir)) dir.create(outdir)
# 2. Plot alignment statistics ----
message("\n[+] Plot alignment statistics \n")
pstats <- plot_alignment_stats(stats)

outfile <- paste0(outdir, "/alignment_statistics.pdf")
message(" -- outfile: ", outfile)
pdf(file = outfile, paper = 'a4r', h = unit(3,'cm'),w = unit(8,'cm'), useDingbats = F)
do.call(gridExtra::grid.arrange, c(pstats, list(nrow=1,ncol=2)))
dev.off()

# 3. Plot gene biotype statistics ----
message("\n[+] Plot gene biotype statistics \n")
pbiotypes <- plot_biotype_stats(stats, colby = colby, biotype = biotype)

for(i in names(pbiotypes)) {
  outfile <- paste0(outdir, "/biotype_stats_",i,".pdf")
  message(" -- outfile: ", outfile)
  pdf(file = outfile, paper = 'a4', h = unit(3,'cm'),w = unit(3,'cm'), useDingbats = F)
  print(pbiotypes[[i]])
  dev.off()
}

message("\n[+] All done. \n")
