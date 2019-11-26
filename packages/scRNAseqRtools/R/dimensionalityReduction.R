# t-distributed Stochastic Neighbor Embedding ----
getTSNE <- function(x
                    , groupBy    = NULL
                    , pal        = NULL
                    , point_size = 1.5
                    , marker     = NULL
                    , ...)
{
  tsne <- Rtsne::Rtsne(x, ...)
  tsneplot <- as.data.frame(tsne$Y)
  rownames(tsneplot) <- rownames(x)
  
  if( !is.null(groupBy) ) {
    if(!is.list(groupBy)) {
      tsneplot$sample <- groupBy[rownames(tsneplot)]
    } else {
      tsneplot$sample     <- groupBy[[1]][rownames(tsneplot)]
      tsneplot$additional <- groupBy[[2]][rownames(tsneplot)]
    }
    
  } else {
    tsneplot$sample <- rownames(tsneplot)
  }
  
  if(is.null(marker)) {
    p0 <- ggplot(tsneplot, aes(x=V1, y=V2, col=sample)) +
      geom_point(size=point_size) + scale_color_manual(values=pal)
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
      geom_point(size=point_size) + facet_wrap(~marker, ncol = 3) +
      scale_color_gradientn(colours = RColorBrewer::brewer.pal(4,"Reds")
                            , guide = guide_colourbar(barheight = 0.6, title.vjust = 1 )
                            , values = c(0.01,0.3,0.6,1)
                            # , limits = c(0,5)
      )
  }
  
  p <- p0 + xlab("t-SNE 1") + ylab("t-SNE 2") +
    theme_bw() + my_theme +
    theme(panel.grid = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1
          , strip.background = element_blank()
          , strip.text = element_text(size = 8, face = "bold"))
  
  
  return(list("tsne" = tsne, "plot" = p))
  
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
