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