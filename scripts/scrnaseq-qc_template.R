# 0. Resources ====
analysis <- "scRNA-seq"

#' ---
#' title: `r eval(analysis)`
#' ---
# command line args ---
args <- commandArgs(trailingOnly = T)

# directories ---
rootdir <- as.character(args[1])
DATADIR   <- paste0(rootdir,"/input_data")
wd <- rootdir

# mapping_stats ---
path_mapping_stats <- paste0(DATADIR,"/multiqc_data")
path_metadata <- paste0(DATADIR,"/metadata.txt")
path_counts <- paste0(DATADIR,"/raw_counts.csv.gz")

# util data ---
# gene_info <- "/Users/andrealauria/epigen/data/datasets/annotations/hg38/gencode.v32.gene.info.tsv"
gene_info <- paste0(DATADIR,"/gencode.v32.gene.info.tsv")


# output dirs ---
FIGDIR  <- paste0(wd,"/Figures/scRNAseq/", analysis)
RESDIR  <- paste0(wd,"/Results/scRNAseq/", analysis)

suppressWarnings(suppressMessages(require(RNAseqRtools)))
suppressWarnings(suppressMessages(require(scRNAseqRtools)))

if(!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = T)
if(!dir.exists(RESDIR)) dir.create(RESDIR, recursive = T)

# helper functions 
get_palette_features <- function(metadata_features)
{
  pals <- ls('package:ggsci', pattern = 'pal')[1:length(metadata_features)]
  palette_features <- lapply(pals, function(p) get(p)()(9))
  names(palette_features) <- metadata_features
  return(palette_features)
}

#' ### 1. Alignment statistics
#' Visualize alignment results
# 1. Alignment statistics ----

metadata <- read.delim(path_metadata, stringsAsFactors = F)
# mapping stats ---
mapping_stats <- read_stats(statdir = path_mapping_stats
                            , metadata = path_metadata
                            , type = 'hisat')

metadata_features <- c('instrument_model','cell.type.ch1')
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

#' ### 2. Biotype statistics 
#' Visualize alignment statistics per transcript categories
# 2. Gene biotype stats ----
biotypes_stats <- calc_gene_biotype_stats(count_matrix = path_counts
                                         , gene_info = gene_info
                                         , metadata = path_metadata
                                         , type = 'hisat'
                                         , coverage_ngenes_th = 10)
#+ fig.width=8, fig.height=8
for(feature in metadata_features) {
  p_biotypes_stats <- plot_all_stats_v2(gene_stats = biotypes_stats
                                        , save_plots = T
                                        , col_by = feature
                                        , th_assigned_reads = 100000
                                        , th_detected_genes = 2000
                                        , th_protein_coding = 75
                                        , th_mitochondrial = 25
                                        , pal = palette_features[[feature]]
                                        , outdir = FIGDIR, w=unit(4,'cm'),h=unit(4,'cm')
                                        , guide_legend_nrow = ceiling(length(unique(biotypes_stats[,feature]))/4))

  
  print(ggpubr::ggarrange(plotlist = p_biotypes_stats
                    , ncol=2, nrow=2
                    , common.legend = TRUE
                    , legend="bottom"))
}

#' ### 3. QC selection
#' Select cells for downstream analysis
# 3. QC selected cells ----
assigned_reads_th <- 60000
ngenes_th <- 2000
perc_mitochondrial_th <- 0

qc_cells <- subset(biotypes_stats, 
                   assigned_reads >= assigned_reads_th &
                     ngenes >= ngenes_th &
                     perc_mitochondrial >= perc_mitochondrial_th)

outfile <- paste0(RESDIR, "/qc_selected_cells.rds")
saveRDS(qc_cells, file = outfile)

outfile <- paste0(RESDIR, "/qc_selected_cells.txt")
write.table(qc_cells, file = outfile, row.names = F, quote = F, sep = "\t")

#' Explore QC results
biotypes_stats$qc_pass <- F
biotypes_stats$qc_pass[match(biotypes_stats$Sample, qc_cells$Sample)] <- T

DT::datatable(biotypes_stats, rownames = F, filter = 'top')