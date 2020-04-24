# Pseudotime Analysis ----
plot_pseudotemporal_ordering <- function(cds, ffeature, cfeature, pal, ...)
{
  #source("theme_setting.R")
  
  if(missing(pal)) pal <- ggsci::pal_aaas()(10)
  
  if(missing(ffeature)) {
    ffeature <- ""
  } else if(missing(cfeature)) {
    cfeature <- ffeature
  } else if(missing(cfeature) & missing(ffeature)) {
    stop(message("[!] Please, provide color/facet feature."))
  }
  
  if(length(ffeature)!=1 || length(cfeature)!=1) {
    
    if(length(cfeature)>1) {
      mergefeature <- cfeature
      cfeature <- "sample"
    } else if(length(ffeature)>1) {
      mergefeature <- ffeature
      ffeature <- "sample"
    }
    
    nm <- apply(phenoData(cds)@data[,mergefeature], 1, paste, collapse = "_")
    phenoData(cds)$sample <- gsub("_","/",nm)
    
    p0 <- plot_cell_trajectory(cds, color_by = cfeature, ...) +
      facet_wrap(paste0("~",ffeature), nrow = 1)
    
  } else {
    p0 <- plot_cell_trajectory(cds, color_by = cfeature, ...)
  }
  
  if(is.numeric(phenoData(cds)@data[,cfeature])) {
    
    p0 <- p0 + scale_color_gradientn(colours = RColorBrewer::brewer.pal(4,"Purples")
                                     , guide = guide_colourbar(barheight = 0.6, title.vjust = 1 )
                                     # , values = c(0.01,0.3,0.6,1)
                                     # , limits = c(-1,1)
    )
  } else {
    p0 <- p0 + scale_color_manual(values = pal)
  }
  pct <- p0 + theme_bw() + my_theme +
    theme(panel.grid = element_blank()
          , strip.background = element_blank()
          , strip.text = element_text(size = 8, face = "bold"))
  
  return(pct)
}

plot_cells_in_pseudotime <- function(cds, structure, pal,point_size=1)
{
  #source("theme_setting.R")
  
  if(length(structure)>1) {
    nm <- apply(pData(cds)[,structure], 1, paste, collapse = "_")
    pData(cds)$sample <- gsub("_","/",nm)
    structure   <- "sample"
  }
  
  ptime  <- pData(cds)[,c("Pseudotime",structure)]
  pptime <- ggplot(ptime, aes_string(x="Pseudotime", y=structure, col=structure)) + geom_point(size=point_size) +
    theme_bw() + my_theme + scale_color_manual(values = pal)
  
  return(pptime)
}

analyze_beam_cluster_profile <- function(pbeam)
{
  #source("theme_setting.R")
  x <- pbeam$heatmap_matrix
  colnames(x) <- c(-99:100)
  beam_clustering <- pbeam$annotation_row
  cp <- lapply(unique(beam_clustering$Cluster),
               function(i)
               {
                 y <- x[rownames(beam_clustering)[beam_clustering$Cluster==i],]
                 z <- apply(y, 2, summary)
                 z <- as.data.frame(t(z))
                 z$time <- as.numeric(rownames(z))
                 colnames(z)[c(2,5)] <- c("qu1","qu3")
                 z$cluster <- i
                 return(list("z"=z,"y"=y))
               } )
  cg <- lapply(cp, function(i) rownames(i$y))
  names(cg) <- paste0("cluster_",unique(beam_clustering$Cluster))
  cp <- lapply(cp, "[[", "z")
  cp <- do.call(rbind.data.frame, cp)
  
  p <- ggplot(cp, aes(x=time, y=Mean, col=cluster)) + geom_line() +
    facet_wrap(~cluster, ncol = 1) + xlab("Pseudotime") + geom_vline(xintercept = 0, lwd=0.25, linetype="dashed") +
    ylab("Smoothed expression") + scale_x_continuous(breaks = 0) +
    geom_ribbon(aes(ymin=qu1, ymax=qu3, fill=cluster), alpha = 0.2) +
    theme_bw() + my_theme + theme(strip.background = element_blank(), strip.text = element_blank())
  
  return(list("p"=p,"gene_clusters"=cg))
}