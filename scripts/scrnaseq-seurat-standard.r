# # # # # # # # # # # # # # # # # # # # # # # #
#                                             #
#  Single cell RNA-sequencing data analysis   #
#                                             #
# # # # # # # # # # # # # # # # # # # # # # # #
#
# scrnaseq-seurat-standard.r
#
# Seurat standard analysis worflow:
#   - QC/filtering
#   - Normalization
#   - Feature selection/scaling
#   - Dimensionality reduction
#   - Clustering
#   - Cluster markers
#   - Visualization
#
# 0. Resources, variables and parameters ----
#
# Mandatory variables to be assigned before running the script:
#   - path_seurat_obj <- ASSIGN_PATH_SEURAT_OBJ
#     path to seurat object with raw counts and metadata (e.g. obtained from scraseq-get_scrna_obj.r script)
#   - path_results    <- ASSIGN_PATH_RESULTS
#     path to results folder
#   - metadata_features <- ASSIGN_METADATA_FEATURES
#     character vector of features in metadata for grouping and visualization
# 
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(scRNAseqRtools)))

# configuration ---
CORES = 12
# Set OpenMP threads
set_openmp_threads = TRUE
if(set_openmp_threads) {
  Sys.setenv("OMP_NUM_THREADS" = CORES)
}

# paths ---
path_seurat_obj <- ASSIGN_PATH_SEURAT_OBJ
path_results    <- ASSIGN_PATH_RESULTS

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# global parameters ---
metadata_features <- ASSIGN_METADATA_FEATURES # used for grouping in plots, character vector

# QC
nFeature_RNA.th = 2000
nCount_RNA.th = 50000
percent.mt.th = 100

# Feature selection and scaling
variable_features_nfeatures <- 2000
variable_features_method <- "vst"
scale_data_features <- NULL
scale_data_vars_to_regress <- NULL

# Dimensionality reduction
max_pcs = 100
npcs = 10
npcs_method = "elbow" # elbow | jackstraw
npc_loadings = 4
npcs_jackstraw = 20
perplexity = 10
tsne_reduction = "pca"
umap_reduction = "pca"

# Clustering
clustering_reduction = "pca"
clsutering_resolution = .8 # seurat default

# Cluster markers 
find_all_markers_only_pos = TRUE
find_all_markers_min_pct  = 0.1
find_all_markers_logfc_th = 0.25
find_all_markers_test_use = "wilcox"
markers_do_heatmap = FALSE
markers_do_vln = TRUE
markers_heatmap_n_top = 20
markers_vln_n_top = 20
markers_vln_ncol = 3
markers_reduction_n_top = 5

# Visualization
reductions = c("pca","tsne","umap")
extended_plots = T # make additional plots for variable features, PCA component scores, markers
reductions_markers = "umap"
reductions_technical_features = "umap"

# read object ---
seurat_obj <- readRDS(path_seurat_obj)

# palette ---
get_palette_features <- function(metadata_features, pals = NULL)
{
  if(is.null(pals)) {
    # set defaults from ggsci
    pals <- ls('package:ggsci', pattern = 'pal')[1:length(metadata_features)] 
  }
  palette_features <- lapply(pals, function(p) {
    # fill missing items fromm RColorBrewer 
    rcb_pal <- sample(rownames(subset(RColorBrewer::brewer.pal.info,category=="div")),1)
    c(get(p)()(9),RColorBrewer::brewer.pal(9,rcb_pal))
  })
  names(palette_features) <- metadata_features
  return(palette_features)
}

set_feature_palettes = T

if(set_feature_palettes) {
  pal_default <- get_palette_features(metadata_features, pals = NULL)
} else {
  pal_default <- ggsci::pal_d3()(10)  
}
# alternatively, set pal_default to manually defined palettes
# pal_default <- ...

rename_seurat_obj = TRUE # rename obj when saving if exists
s33d = 1991
# 1. QC/filtering ----
# clean object
# remove zero count genes in all cells
seurat_obj <- seurat_obj[-which(rowSums(seurat_obj)==0),]

# remove empty/minibulk wells
cell_rm_idx <- grep("EMPTY|(MINI)?BULK",colnames(seurat_obj),ignore.case = T)
if(length(cell_rm_idx)!=0) seurat_obj <- seurat_obj[,-cell_rm_idx]

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^(MT|mt)-")

# qc plots
qc_vlnplot <- list()
for(feature in metadata_features) {
  qc_vlnplot[[feature]] <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = feature) &
    theme_bw() + my_theme + theme(legend.key.size=unit(4,'mm'), axis.text.x = element_text(angle=45,hjust=1,vjust=1)) &
    scale_fill_manual(values = pal_default[[feature]])
  
  outfile <- paste0(path_results,"/qc_vlnplot.",feature,".pdf")
  message(" -- output file: ", outfile)
  pdf(file = outfile, paper = "a4r", useDingbats = F, w=unit(10,'cm'), height = unit(4,'cm'))
  print(qc_vlnplot[[feature]])
  dev.off()
}

# filter
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= nFeature_RNA.th & nCount_RNA >= nCount_RNA.th & percent.mt < percent.mt.th)

# 2. Normalization ----
seurat_obj <- NormalizeData(seurat_obj)

# 3. Feature selection and scaling ----
seurat_obj <- FindVariableFeatures(object = seurat_obj
                                   , selection.method = variable_features_method
                                   , nfeatures = variable_features_nfeatures)

if(extended_plots) {
  pvarfeatures <- VariableFeaturePlot(seurat_obj, pt.size = 0.25) + theme_bw() + my_theme
  pvarfeatures <- LabelPoints(plot = pvarfeatures
                              , points = head(VariableFeatures(seurat_obj), 10)
                              , repel = TRUE, xnudge = 0, ynudge = 0, size=2)
  
  outfile <- paste0(path_results,"/variable_features.pdf")
  message(" -- output file: ", outfile)
  pdf(file = outfile, paper = "a4r", useDingbats = F, w=unit(3.5,'cm'), height = unit(3.5,'cm'))
  print(pvarfeatures)
  dev.off()
}


seurat_obj <- ScaleData(seurat_obj
                        , features = scale_data_features
                        , vars.to.regress = scale_data_vars_to_regress)

# 4. Dimensionality reduction ----
# PCA
seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = max_pcs)

if("jackstraw"%in%npcs_method) {
  seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
  seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:npcs_jackstraw)
  pjackstraw <- JackStrawPlot(seurat_obj, dims = 1:npcs_jackstraw) + 
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'))
  
  outfile <- paste0(path_results, "/PCs_dimensionality.jackstraw.pdf")
  pdf(file = outfile, paper = "a4", w = unit(3.5, 'cm'), h = unit(4, 'cm'), useDingbats = F)
  print(pjackstraw)
  dev.off()
}
if("elbow"%in%npcs_method) {
  pelbow <- ElbowPlot(seurat_obj) + 
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'))
  
  outfile <- paste0(path_results, "/PCs_dimensionality.elbow.pdf")
  pdf(file = outfile, paper = "a4", w = unit(3.5, 'cm'), h = unit(4, 'cm'), useDingbats = F)
  print(pelbow)
  dev.off()
}

ppca_loadings <- VizDimLoadings(seurat_obj, dims = 1:npc_loadings, reduction = "pca") & theme_bw() + my_theme
outfile <- paste0(path_results, "/PCs_loadings.pdf")
message(" -- writing to: ", outfile)
pdf(file = outfile, paper = "a4", w = unit(10, 'cm'), h = unit(10, 'cm'), useDingbats = F)
print(ppca_loadings)
dev.off()

if(extended_plots) {
  outfile <- paste0(path_results, "/PCs_loadings.heatmap.pdf")
  message(" -- writing to: ", outfile)
  pdf(file = outfile, paper = "a4", w = unit(10, 'cm'), h = unit(10, 'cm'), useDingbats = F)
  DimHeatmap(seurat_obj, dims = 1:npc_loadings, cells = 500, balanced = TRUE, ncol = 2)
  dev.off()
}
# tSNE
set.seed(s33d)
seurat_obj <- Seurat::RunTSNE(seurat_obj
                              , reduction = tsne_reduction
                              , dims = 1:npcs
                              , perplexity=perplexity)
# UMAP
set.seed(s33d)
seurat_obj <- Seurat::RunUMAP(seurat_obj
                              , reduction = umap_reduction
                              , dims = 1:npcs)
# 5. Clustering ----
seurat_obj <- FindNeighbors(seurat_obj
                            , reduction = clustering_reduction
                            , dims = 1:npcs)
seurat_obj <- FindClusters(seurat_obj
                           , resolution = clsutering_resolution)

# 6. Find cluster markers ----
markers <- FindAllMarkers(seurat_obj
                          , only.pos = find_all_markers_only_pos
                          , min.pct = find_all_markers_min_pct
                          , logfc.threshold = find_all_markers_logfc_th
                          , test.use = find_all_markers_test_use)

outfile <- paste0(path_results,"/seurat_clusters_markers.txt.gz")
write.table(markers, file = gzfile(outfile), row.names = F, col.names = T, sep = "\t", quote = F)

# 7. Visualization ----
dimplot_wrapper <- function(seurat_obj
                            , metadata_features
                            , pal = pal_default
                            , reduction = "pca") 
{
  pdimred_list <- list()
  if(!is.list(pal)) {
    pal_dimred <- list("seurat_clusters"=pal)
  } else {
    pal_dimred <- pal
  }
  pdimred_list[["seurat_clusters"]] <- DimPlot(seurat_obj, reduction = reduction, group.by = "seurat_clusters", pt.size = .5) +
    theme_bw() + my_theme + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'), aspect.ratio = 1) +
    xlab(paste0(gsub("PCA","PC",toupper(reduction))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reduction))," 2")) + ggtitle("seurat_clusters") + scale_color_manual(values = pal_dimred[["seurat_clusters"]])
  
  for(feature in metadata_features) {
    if(!is.list(pal)) pal_dimred[[feature]] <- pal
    pdimred_list[[feature]] <- DimPlot(seurat_obj, reduction = reduction, group.by = feature, pt.size = .5) +
      theme_bw() + my_theme + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'), aspect.ratio = 1) +
      xlab(paste0(gsub("PCA","PC",toupper(reduction))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reduction))," 2")) + ggtitle(feature) + scale_color_manual(values = pal_dimred[[feature]])
  }
  
  pdimred <- ggpubr::ggarrange(plotlist = pdimred_list, align = "hv",nrow=1)
  pdimred <- ggpubr::annotate_figure(pdimred, top=text_grob(paste0("Dimensionality reduction - ",toupper(reduction)), hjust=0.5, size=8))
  
  return(pdimred)
}

# palettes
if(set_feature_palettes) {
  pal_cluster <- get_palette_features(metadata_features = "seurat_clusters"
                       ,pals = ls('package:ggsci', pattern = 'pal')[length(metadata_features)+1])
  pal_default <- c(pal_default, pal_cluster)
}

pdimred <- list()
for(reduction in reductions) {
  pdimred[[reduction]] <- dimplot_wrapper(seurat_obj, metadata_features, reduction = reduction, pal = pal_default)
  
  outfile <- paste0(path_results, "/dimplot.",reduction,".pdf")
  message(" -- saving to: ",outfile)
  pdf(file = outfile, paper = "a4r", w = unit(10, 'cm'), h = unit(6, 'cm'), useDingbats = F)
  print(pdimred[[reduction]])
  dev.off()  
}

# markers
if(markers_do_vln) {
  
  top_n_genes <- lapply(split(markers, markers$cluster), function(x) head(x, n = markers_vln_n_top))
  top_n_genes <- do.call(rbind.data.frame, top_n_genes)

  if(length(top_n_genes$gene)>markers_vln_ncol*5) {
    top_n_genes_list <- split(top_n_genes$gene, ceiling(seq_along(top_n_genes$gene)/as.integer(markers_vln_ncol*5)))
  } else {
    top_n_genes_list <- list(top_n_genes$gene)
  }
  pmarkers_vln <- list()
  for(i in 1:length(top_n_genes_list)) {
    pmarkers_vln[[i]] <- VlnPlot(seurat_obj, features = top_n_genes_list[[i]]
                                 , group.by = "seurat_clusters"
                                 , cols = pal_default[["seurat_clusters"]]
                                 , pt.size = 0.05,ncol = markers_vln_ncol) &
      theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'))
    
  }
  
  outfile <- paste0(path_results, "/seurat_clusters_markers.vln.top_",markers_vln_n_top,".pdf")
  message(" -- writing to: ", outfile)
  pdf(file = outfile, paper = "a4", w = unit(20, 'cm'), h = unit(20, 'cm'), useDingbats = F)
  print(pmarkers_vln)
  dev.off()
 
}
if(markers_do_heatmap) {
  top_n_genes <- lapply(split(markers, markers$cluster), function(x) head(x, n = markers_heatmap_n_top))
  top_n_genes <- do.call(rbind.data.frame, top_n_genes)
  
  outfile <- paste0(path_results, "/seurat_clusters_markers.heatmap.top_",markers_heatmap_n_top,".pdf")
  pdf(file = outfile, paper = "a4", w = unit(8, 'cm'), h = unit(10, 'cm'), useDingbats = F)
  DoHeatmap(seurat_obj, features = top_n_genes$gene) + NoLegend()
  dev.off()
}

if(extended_plots) {
  # technical feature plots
  pf1 <- FeaturePlot(seurat_obj, features = "nFeature_RNA",reduction = reductions_technical_features, pt.size = .5) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm')) +
    xlab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 2")) 
  pf2 <- FeaturePlot(seurat_obj, features = "nCount_RNA",reduction = reductions_technical_features, pt.size = .5) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm')) +
    xlab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 2")) 
  pf3 <- FeaturePlot(seurat_obj, features = "percent.mt",reduction = reductions_technical_features, pt.size = .5) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm')) +
    xlab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 2")) 
  
  pdimred_technical_features <- ggpubr::ggarrange(pf1,pf2,pf3,align = "hv", nrow = 1, ncol = 3)
  outfile <- paste0(path_results, "/dimplot.",reductions_technical_features,"_technical_features.pdf")
  pdf(file = outfile, paper = "a4r", w = unit(10, 'cm'), h = unit(4, 'cm'), useDingbats = F)
  print(pdimred_technical_features)
  dev.off()
  
  top_n_genes <- lapply(split(markers, markers$cluster), function(x) head(x, n = markers_reduction_n_top))
  top_n_genes <- do.call(rbind.data.frame, top_n_genes)
  
  pdimred_markers <- list()
  for(reduction in reductions_markers) {
    pdimred_markers[[reduction]] <- FeaturePlot(seurat_obj, features = top_n_genes$gene, reduction = reduction, pt.size = .1, ncol = 5) &
      theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'))
    
    outfile <- paste0(path_results,"/seurat_clusters_markers.dimred_",reduction,".pdf")
    message(" -- output file: ", outfile)
    pdf(file = outfile, paper = "a4", useDingbats = F, w=unit(30,'cm'), height = unit(30,'cm'))
    print(pdimred_markers[[reduction]])
    dev.off()
  }
}
# 8. Save object ----
outfile <- paste0(path_results,"/seurat_obj.rds")
if(file.exists(outfile) && rename_seurat_obj) {
  message("[!] Object exists, renaming")
  outfile <- gsub("\\.rds$",".proc.rds",outfile)
}
message(" -- saving seurat object to: ", outfile)
saveRDS(seurat_obj, file = outfile)
