#!/usr/bin/env Rscript

# # # # # # # # # # # # # # # # # # # # # # # # # 
# single cell RNA-sequencing - quality controls #
# # # # # # # # # # # # # # # # # # # # # # # # #

# 0. Resources ----
# Command line/Env variables ---
suppressWarnings(suppressMessages(require(docopt)))
'Usage:
   scrnaseq-qc.r [-q <qcmatrix> -c <counts> -m <metadata> -g <gencode> -p <pipeline> -b <biotype> -f <features> -a <analysis> -o <outdir>]

Options:
   -q, --qcmatrix Path to QC matrix or multiQC output directory.
   -c, --counts Path to raw counts. 
   -m, --metadata Path to experiment metadata.
   -g, --gencode Path to GENCODE gene/transcript biotype annotation table.
   -p, --pipeline Pipeline used for data processing (hisat/salmon/star) [default: hisat].
   -b, --biotype Gene biotype (comma separated if > 1) for diagnostic plots. 
                 Protein coding, lncRNA, rRNA and mitochondrial genes 
                 are reported by default.
   -f, --features Color cells by feature (comma separated if > 1) in metadata.
   -a, --analysis Project title [default: scRNA-seq].
   -o, --outdir Output directory. [default: .]
' -> doc

opts <- docopt(doc)
required_args <- opts[c(1:4,7)]
if(any(sapply(required_args, is.null))) {
  missing_idx <- sapply(required_args, is.null)
  missing_args <- gsub("--"," ",names(required_args)[missing_idx])
  message("\n[!] Missing required arguments: ",missing_args,"\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

# mapping_stats ---
path_mapping_stats <- as.character(opts$qcmatrix)
path_counts   <- as.character(opts$counts)
path_metadata <- as.character(opts$metadata)

# util data ---
biotype  <- opts$biotype
gencode  <- as.character(opts$gencode)
pipeline <- as.character(opts$pipeline)
metadata_features <- as.character(opts$features)
outdir   <- as.character(opts$outdir)
analysis <- as.character(opts$analysis)

# directories ---
FIGDIR <- paste0(outdir,"/Figures/", analysis)
RESDIR <- paste0(outdir,"/Results/", analysis)

if(!dir.exists(outdir)) dir.create(outdir)
if(!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = T)
if(!dir.exists(RESDIR)) dir.create(RESDIR, recursive = T)

# requirements ---
suppressWarnings(suppressMessages(require(scRNAseqRtools)))
# helper functions 
get_comma_arglist <- function(argument)
{
  if(!is.null(argument) && grepl("\\,",argument)) {
    argument <- unlist(strsplit(argument,"\\,"))
  }
  return(argument)
}
get_palette_features <- function(metadata_features)
{
  pals <- ls('package:ggsci', pattern = 'pal')[1:length(metadata_features)]
  palette_features <- lapply(pals, function(p) get(p)()(9))
  names(palette_features) <- metadata_features
  return(palette_features)
}

#' ---
#' title: `r eval(analysis)`
#' ---
#' ### 1. Alignment statistics
#' Visualize reads alignment results, summarize by sample features.
# 1. Alignment statistics ----
metadata_features <- get_comma_arglist(metadata_features)
biotype <- get_comma_arglist(biotype)

# mapping stats ---
mapping_stats <- read_stats(statdir = path_mapping_stats
                            , metadata = path_metadata
                            , pipeline = pipeline)

palette_features  <- get_palette_features(metadata_features)

#+ fig.width=8, fig.height=8
for(feature in metadata_features) {
  
  p_map_stats_feature <- plot_mapping_stats(mapping_stats
                                            , color_by = feature
                                            , pal = palette_features[[feature]]
                                            , p_boxplot = T)
  print(p_map_stats_feature)
  outfile <- paste0(FIGDIR,"/Mapping_rate_colby_",feature,".pdf")
  message(" -- output file: ", outfile)
  pdf(file = outfile, paper = "a4r", useDingbats = F, w=unit(6,'cm'), height = unit(6,'cm'))
  print(p_map_stats_feature)
  dev.off()
  
}

#' Overall alignment distribution
p_overall_stat <- plot_overall_alignment_stats(mapping_stats = mapping_stats)
#+ fig.width=5, fig.height=4
p_overall_stat
outfile <- paste0(FIGDIR,"/overall_mapping_rate.pdf")
pdf(file = outfile, paper = "a4r", useDingbats = F, w=unit(3.5,'cm'), height = unit(2.5,'cm'))
p_overall_stat
dev.off()

#' ### 2. Gene biotype statistics 
#' Genes/transcripts QC alignment results. Visualize alignment statistics per transcript categories and sample features.
# 2. Gene biotype stats ----
biotypes_stats <- calc_gene_biotype_stats(count_matrix = path_counts
                                         , gene_info = gencode
                                         , biotypes = biotype
                                         , metadata = path_metadata
                                         , pipeline = pipeline
                                         , coverage_ngenes_th = 10)
print(head(biotypes_stats))
#+ fig.width=9, fig.height=20
for(feature in metadata_features) {
  p_biotypes_stats <- plot_all_stats_v3(gene_stats = biotypes_stats
                                        , biotypes = biotype
                                        , save_plots = T
                                        , col_by = feature
                                        , th_assigned_reads = 100000
                                        , th_detected_genes = 2000
                                        , th_mitochondrial = 25
                                        , pal = palette_features[[feature]]
                                        , outdir = FIGDIR, w=unit(7,'cm'),h=unit(4,'cm')
                                        , guide_legend_nrow = ceiling(length(unique(biotypes_stats[,feature]))/4))

  
  print(ggpubr::ggarrange(plotlist = p_biotypes_stats
                    , ncol=1
                    , common.legend = TRUE
                    , legend="none"))
}

#' ### 3. Explore QC results
#' Global alignment statistics
DT::datatable(mapping_stats, rownames = F, filter = 'top')

outfile <- paste0(RESDIR, "/mapping_stats.rds")
saveRDS(mapping_stats, file = outfile)

outfile <- paste0(RESDIR, "/mapping_stats.txt")
write.table(mapping_stats, file = outfile, row.names = F, quote = F, sep = "\t")

#' Gene biotype statistics
DT::datatable(biotypes_stats, rownames = F, filter = 'top')

outfile <- paste0(RESDIR, "/biotypes_stats.rds")
saveRDS(biotypes_stats, file = outfile)

outfile <- paste0(RESDIR, "/biotypes_stats.txt")
write.table(biotypes_stats, file = outfile, row.names = F, quote = F, sep = "\t")
