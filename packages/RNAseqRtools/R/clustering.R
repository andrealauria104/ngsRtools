# Tools for RNA-seq data analysis


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
