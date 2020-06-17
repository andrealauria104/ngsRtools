# t-distributed Stochastic Neighbor Embedding ----
getTSNE <- function(x
                    , groupBy = NULL
                    , experimental_info = NULL
                    , col_by            = NULL
                    , shape_by          = NULL
                    , pal               = NULL
                    , point_size        = 1.5
                    , marker            = NULL
                    , facet_ncol        = 3 
                    , scale_col_limits  = c(-1.5,1.5)
                    , ... )
{
  tsne <- Rtsne::Rtsne(x, ...)
  tsneplot <- as.data.frame(tsne$Y)
  rownames(tsneplot) <- rownames(x)
  
  if(!is.null(groupBy)) experimental_info <- groupBy # for compatibility with previous version
  
  if( !is.null(experimental_info) ) {
    if(is.character(experimental_info)) {
      # for compatibility with previous version
      experimental_info <- reshape2::melt(experimental_info, value.name="sample")
      if(!is.null(col_by)) {
        colnames(experimental_info)[1] <- col_by 
      }
    }
    if(is.null(col_by)) {
      col_by <- colnames(experimental_info)[1]
    }
    tsneplot[[col_by]] <- experimental_info[rownames(tsneplot), col_by]
    if(!is.null(shape_by)) {
      tsneplot[[shape_by]] <- experimental_info[rownames(tsneplot), shape_by]
    }
  } else {
    if(is.null(col_by)) col_by <- "sample"
    tsneplot[[col_by]] <- rownames(tsneplot) 
  }
  
  if(is.null(pal)) {
    pal <- ggsci::pal_aaas()
    pal <- pal(length(unique(tsneplot[,col_by])))
  }
  
  if(is.null(marker)) {
    print(head(tsneplot))
    p0 <- ggplot(tsneplot, aes_string(x="V1", y="V2", col=col_by, shape=shape_by)) +
      geom_point(size=point_size) + scale_color_manual(values=pal)+
      guides(col = guide_legend(ncol=1))
    
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
      geom_point(size=point_size) + facet_wrap(~marker, ncol = facet_ncol) +
      scale_colour_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B"
                             , midpoint = 0
                             , breaks = seq(-2,2,1)
                             , limits = scale_col_limits
                             , oob=squish) +
      guides(col = guide_colourbar(barwidth = 0.8))
      
  }
  
  p <- p0 + xlab("t-SNE 1") + ylab("t-SNE 2") +
    theme_bw() + my_theme +
    theme(panel.grid = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1
          , strip.background = element_blank()
          , strip.text = element_text(size = 8, face = "bold")
          , legend.key.size =  unit(0.5,'cm')
          , legend.position = 'right')
  
  
  return(list("tsne" = tsne, "plot" = p))
  
}

prepare_tsneplot_genes <- function(sce, genes
                                   , cumulative_expression = F
                                   , assay.var = "normcounts")
{
  if(!cumulative_expression || length(genes)==1) {
    m_genes <- t(apply(assay(sce, assay.var)[genes,, drop=F],1, scale))
    colnames(m_genes) <- colnames(sce)
    tsneplot_genes <- reshape2::melt(m_genes
                                     , varnames = c("name","cell")
                                     , value.name = "expression")
  } else {
    geneset <- lapply(genes, function(x) {
      idx <- intersect(x, rownames(sce))
      colSums(assay(sce, assay.var)[idx,])
    })
    geneset <- do.call(rbind, geneset)
    geneset
    m_geneset <- t(apply(geneset, 1, scale))
    colnames(m_geneset) <- colnames(geneset)
    tsneplot_genes <- reshape2::melt(m_geneset
                                     , varnames = c("name","cell")
                                     , value.name = "expression")
  }
  
  
  return(tsneplot_genes)
}

# Dimensionality reduction ----
plot_dimred_cells <- function(cds, pal, ...)
{
  if(missing(pal)) pal <- 'grey'
  p <- plot_cells(cds, ...) +
    scale_color_manual(values = pal) +
    theme_bw() + my_theme + theme(panel.grid = element_blank())
  
  return(p)
}
