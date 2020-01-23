# Gene Ontology ====
getGO <- function(gene_list
                  , ORGANISM = "mmusculus"
                  , ordered_query = F
                  , correction_method='fdr'
                  , hier_filtering ='moderate'
                  , min_set_size = 10
                  , process = T
                  , pvth = 0.05
                  , ncut = NULL
                  , ...)
{
 
  src_filters = c("GO:BP","GO:MF","GO:CC","KEGG",
                 "REAC","TF","MI","CORUM","OMIM", 
                 "HPA","HP")
  
  profileres <- lapply(src_filters, function(src)
    {
      gProfileR::gprofiler( gene_list,
                            organism = ORGANISM,
                            max_p_value = 0.05,
                            ordered_query = ordered_query,
                            significant = F,
                            correction_method = correction_method,
                            min_set_size = min_set_size,
                            hier_filtering = hier_filtering,
                            src_filter = src,
                            ... )}
    )
  names(profileres) <- src_filters
  
  if(process) profileres <- lapply(profileres, processGO, pvth = pvth, ncut = ncut)
  
  return(profileres)
}

processGO <- function(x, pvth = 0.05, ncut = NULL)
{
  x <- subset(x, p.value<=pvth)
  x <- x[,c("term.id","term.name","domain","p.value","intersection")]
  if( any(duplicated(x$term.name))) x <- x[-which(duplicated(x$term.name)),]
  x$term.name <- factor(x$term.name, 
                        levels = x$term.name[order(x$p.value, decreasing = F)])
  if(!is.null(ncut)){
    if(nrow(x)>ncut) x <- x[1:ncut,]
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
    ggplot(go, aes(x=Description, y=-log10(qvalue))) + geom_col(position = "dodge", alpha = 0.8, ...) +
      coord_flip() + theme_bw() + my_theme_2 + theme(axis.title.y=element_blank())
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
    db <- org.Mm.eg.db::org.Mm.eg.db
  } else if(species=='hs') {
    db <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    stop(message("[!] Provide a valid organism (available: mm/hs)."))
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
    ego <- clusterProfiler::simplify(ego)
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
        colnames(x)[grep("term.name|p.value",colnames(x))] <- c("Description","p.adjust")
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
                           , row_names_side = "left"
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
                                  , column_names_gp = gpar(fontsize=term_font_size)
                                  # , row_names_max_width = unit(5, "cm")
                                  # , width = unit(2,'cm')
                                  , row_names_side = row_names_side 
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

process_enrich_results <- function(enrich_res, gene_info
                                   , fdrth = NULL
                                   , dea = NULL
                                   , outfile = NULL
                                   , outplot = NULL)
{
  if(any(grep("enrichResult",class(enrich_res)))) {
    enrich_res <- enrich_res@result
  }
  
  if(any(grep("mm", enrich_res$ID))) {
    enrich_res$geneID <- sapply(enrich_res$geneID, function(x)
    {
      x <- unlist(strsplit(x,"/"))
      nm <- as.character(gene_info$external_gene_name[match(unlist(strsplit(x,"/")), gene_info$entrezgene)])
      nm <- paste0(nm,collapse = "/")
      return(nm)
    })
  }
  
  if(!is.null(deg)) {
    if(any(grep("TopTags",class(deg)))) {
      deg <- deg$table
    }
    nde <- lapply(enrich_res$geneID, function(x)
    {
      x <- unlist(strsplit(x,"/"))
      idx <- match(x, deg$genes)
      nup <- sum(deg[idx,'logFC']>0)
      ndw <- sum(deg[idx,'logFC']<0)
      return(data.frame("N.up" = nup,"N.dw"=ndw))
    })
    enrich_res <- cbind.data.frame(enrich_res, nde)
  }
  if(!is.null(fdrth)) {
    enrich_res <- subset(enrich_res, qvalue <= fdrth)
  }
  if(!is.null(outfile)) {
    saveXLSresEdgeR(enrich_res, outfile = outfile, name = "KEGG enrichment results")
  }
  
  if(!is.null(outplot)) {
    toplot <- enrich_res[,c("Description","N.up","N.dw","qvalue")]
    toplot$descore <- (toplot$N.up-toplot$N.dw)*(-log10(toplot$qvalue))
    toplot <- toplot[order(abs(toplot$descore), decreasing = T),]
    nm <- toplot$Description
    toplot$Description <- factor(toplot$Description, levels = nm)
    toplot$direction <- factor(sign(toplot$descore))
    p <- ggplot(toplot[1:15,], aes(y=descore, x=Description, fill = direction)) + geom_col(alpha=0.7) +
      theme_bw() + my_theme_2 + scale_fill_manual(values = c('darkblue','darkred', 'lightgrey')) +
      coord_flip()
    
    pdf(file = outplot, paper = 'a4', w = unit(4,'cm'), h = unit(4,'cm'), useDingbats = F)
    print(p)
    dev.off()
    
  }
}