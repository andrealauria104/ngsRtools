# Principal Component Analysis ====
getPCA <- function(x, cor.test.loadings=F) 
{
  message("[+] Run PCA ...")
  rna <- t(x)
  rna <- rna[,colSums(rna)>1]
  pca <- prcomp(rna, scale. = T, center = T)
  pca_summary <- summary(pca)$importance
  pca_load <- t(t(pca$rotation)*pca_summary[1,])
  
  pca_data <- list("pca" = pca, "pca_summary" = pca_summary, "pca_load" = pca_load)
  
  if(cor.test.loadings) {
    message("[+] Variable/components correlation testing ...")
    pca_cor.pv <- matrix(data = NA, nrow=ncol(rna), ncol = nrow(pca$x)
                         , dimnames = list(colnames(rna), colnames(pca$x)))
    for(i in 1:ncol(rna)) {
      for(j in 1:nrow(pca$x)) {
        pca_cor.pv[i,j] <- cor.test(rna[,i],pca$x[,j])$p.value
      }
    }
    pca_data$pca_cor.pv <- pca_cor.pv
  }
  
  return(pca_data)
}

plotPCA <- function(x, experimental_info = NULL
                    , labels   = T
                    , col_by   = NULL
                    , shape_by = NULL
                    , pal      = NULL
                    , scree_plot = F
                    , scree_plot_type = "explained_variance"
                    , point_size = 2
                    , max.pcs = NULL
                    , dim_1 = "PC1"
                    , dim_2 = "PC2"
                    , default.title = "RNA-seq PCA"
                    , legend_position = "bottom"
                    , col_legend_nrow = NULL
                    , col_legend_ncol = NULL
                    , shape_legend_nrow = NULL
                    , shape_legend_ncol = NULL)
{
  if(is.list(x) && any(grepl("pca",names(x)))) {
    pca_data <- x
  } else {
    pca_data <- getPCA(x)
  }
  
  if((is.null(max.pcs) || max.pcs > ncol(pca_data$pca$x))) {
    max.pcs <- ncol(pca_data$pca$x)
  }
  
  pca_plot <- as.data.frame(pca_data$pca$x[,paste0("PC",1:max.pcs)])
  
  if( !is.null(experimental_info) ) {
    if(is.character(experimental_info)) {
      # for compatibility with previous version
      experimental_info <- reshape2::melt(experimental_info)
      if(!is.null(col_by)) {
        colnames(experimental_info)[1] <- col_by 
      }
    }
    if(is.null(col_by)) {
      col_by <- colnames(experimental_info)[1]
    }
    pca_plot[[col_by]] <- experimental_info[rownames(pca_plot), col_by]
    if(!is.null(shape_by)) {
      pca_plot[[shape_by]] <- experimental_info[rownames(pca_plot), shape_by]
    }
  } else {
    if(is.null(col_by)) col_by <- "samples"
    pca_plot[[col_by]] <- rownames(pca_plot) 
  }
  
  if(is.null(pal)) {
    pal <- ggsci::pal_d3()
    pal <- pal(length(unique(pca_plot[,col_by])))
  }
  pca_plot$repel_col_by <- pca_plot[,col_by]
  
  p <- ggplot(pca_plot, aes_string(x=dim_1, y=dim_2, col=col_by, shape=shape_by)) + geom_point(size=point_size) +
    xlab(paste0(dim_1," (",round(pca_data$pca_summary[2,grep(paste0(dim_1,"$"), colnames(pca_data$pca_summary))]*100,1),"%)")) +
    ylab(paste0(dim_2," (",round(pca_data$pca_summary[2,grep(paste0(dim_2,"$"), colnames(pca_data$pca_summary))]*100,1),"%)")) + 
    theme_bw() + ggtitle(default.title) + my_theme_2 +
    theme(panel.grid.minor = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1
          , legend.key.size = unit(5,'mm')
          , legend.position = legend_position) +
    scale_color_manual(values=pal)
  
  if(!is.null(col_legend_nrow)) {
    p <- p + guides(col = guide_legend(nrow=col_legend_nrow))
  } 
  if(!is.null(col_legend_ncol)) {
    p <- p + guides(col = guide_legend(ncol=col_legend_ncol)) 
  }
  if(!is.null(shape_legend_nrow)) {
    p <- p + guides(shape = guide_legend(nrow=shape_legend_nrow))
  } 
  if(!is.null(shape_legend_ncol)) {
    p <- p + guides(shape = guide_legend(ncol=shape_legend_ncol)) 
  }
  
  if(labels) {
    p <- p + ggrepel::geom_label_repel(aes(label = rownames(pca_plot), col=repel_col_by),
                                       fontface = 'bold'
                                       # , color = 'black'
                                       , show.legend = F
                                       , size=2,
                                       box.padding = 0.35, point.padding = 0.5,
                                       segment.color = 'grey50') 
  }
  
  if(scree_plot) {
    
    if(scree_plot_type=='explained_variance') {
      spdata      <- reshape2::melt(pca_data$pca_summary[c(2:3),])
      spdata$Var2 <- as.numeric(gsub("PC", "",spdata$Var2))
      spdata      <- spdata[spdata$Var2<=max.pcs,]
      sp <- ggplot(spdata, aes(x=Var2, y=value, group=Var1, col=Var1)) + geom_point(size=1) + geom_line() +
        theme_bw() + my_theme_2 + scale_x_continuous(breaks = spdata$Var2) + 
        geom_hline(yintercept = 0.9, lwd=0.2, linetype="dashed", col = "darkred") +
        xlab("Principal Component") + ylab("Explained Variance") + ggtitle("PCA - Scree Plot") + 
        theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank()) +
        scale_y_continuous(breaks = seq(0,1,0.1),labels = scales::percent, limits = c(0,1)) + scale_color_manual(values=c("black","grey"))
    } else if(scree_plot_type=='standard_deviation') {
      spdata      <- reshape2::melt(pca_data$pca_summary[1,])
      
      spdata$Var2 <- as.numeric(gsub("PC", "",rownames(spdata)))
      spdata      <- spdata[spdata$Var2<=max.pcs,]  
      
      sp <- ggplot(spdata, aes(x=Var2, y=value), col='black') + geom_point(size=1) + geom_line() +
        theme_bw() + my_theme_2 + scale_x_continuous(breaks = spdata$Var2) + 
        geom_hline(yintercept = 0.9, lwd=0.2, linetype="dashed", col = "darkred") +
        xlab("Principal Component") + ylab("Standard Deviation") + ggtitle("PCA - Scree Plot") + 
        theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank())
      
    }
    
    return(list("pca" = p, "screeplot" = sp))
  } else {
    return(p) 
  }
}

assignLoadings <- function(pca_data
                           , n.comp = NULL
                           , cor.test.loadings = F
                           , cor.pv.th = 0.05)
{
  if(is.list(pca_data)) {
    if(!any(grepl("pca_load",names(pca_data$pca)))) {
      pca_load <- t(t(pca_data$pca$rotation)*pca_data$pca_summary[1,])
    } else {
      pca_load <- pca$pca_load
    }
  } else {
    pca_load <- pca_data
  }
  
  assigned <- apply(pca_load, 1, function(x) which.max(abs(x)))
  assigned <- paste0("PC",assigned)
  assigned <- data.frame("variable"    = rownames(pca_load)
                         , "princomp"  = assigned
                         , stringsAsFactors = F)
  
  if(cor.test.loadings) {
    message(" -- Flter genes/loadings by significant correlation")
    message(" -- P-value threshold: ", cor.pv.th)
    pvs <- vector()
    if(any(grepl("pca_cor.pv",names(pca_data)))) {
      for(i in 1:nrow(assigned)) {
        gene <- assigned[i,1]
        pc <- assigned[i,2]
        pvs[gene] <- pca_data$pca_cor.pv[i,pc]
      }
      assigned$pv <- pvs[assigned$variable]
      assigned <- subset(assigned, pv<=cor.pv.th)
    } else {
      warning("[!] Missing p-values from cor.test. Re-run getPCA with cor.test.loadings=TRUE for gene filtering.")
    }
  }
  
  if(!is.null(n.comp)) {
    message(" -- Keeping n = ",n.comp, " Principal Components")
    if( is.numeric(n.comp) )  {
      if( length(n.comp)==1 )  n.comp <- 1:n.comp  
      n.comp <- paste0("PC",n.comp)
    }
    assigned <- subset(assigned, princomp%in%n.comp)
  }
  
  return(assigned)
}

plotLoadings <- function(pca_data, components
                         , point_size = 0.2
                         , assigned = NULL
                         , genes = NULL
                         , genes_info = NULL
                         , genes.jitter.width = 0.1
                         , genes.jitter.height = 0.15)
{
  # draw radius-1 circle for loadings/correlations plot
  circle <- function(center = c(0, 0), npoints = 100) {
    r = 1
    tt = seq(0, 2 * pi, length = npoints)
    xx = center[1] + r * cos(tt)
    yy = center[1] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  corcir <- circle(c(0, 0), npoints = 100)
  
  toplot <- as.data.frame(pca_data$pca_load[,components])
  
  if(!is.null(assigned)) {
    idx <- intersect(assigned$variable, rownames(toplot))
    toplot <- toplot[idx, ]
    aidx <- match(assigned$variable, rownames(toplot))
    toplot$assigned <- assigned[aidx, 'princomp']
    nidx <- which(!toplot$assigned%in%components)
    toplot$assigned[nidx] <- "others"
    
    p0 <- ggplot(toplot, aes_string(x=components[1], y=components[2], col="assigned"))
  } else {
    p0 <- ggplot(toplot, aes_string(x=components[1], y=components[2]))
  }
  
  p <- p0 + geom_point(size=point_size, alpha=0.8) + 
    geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") +
    geom_hline(yintercept = 0, lwd=0.2, linetype="dashed", col = "black") +
    geom_vline(xintercept = 0, lwd=0.2, linetype="dashed", col = "black") +
    ggtitle("PCA Loadings Plot") + scale_color_manual(values = c("#404040","#CC0000","orange")) +
    theme_bw() + my_theme + theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank()) 
  
  if(!is.null(genes)) {
    if(is.null(genes_info)) genes_info <- rep("Gene loading vector",length(genes))
    arrows <- data.frame(x1 = rep(0,length(genes)), y1 = rep(0,length(genes))
                         , x2 = pca_data$pca_load[genes,components[1]]
                         , y2 = pca_data$pca_load[genes,components[2]]
                         , gene_name = genes
                         , genes_info = genes_info)
    
    p.genes <- ggplot(toplot, aes_string(x=components[1], y=components[2])) + 
      geom_text(data=arrows, aes(x = x2, y = y2, label=gene_name), size = 2, position=position_jitter(width=genes.jitter.width,height=genes.jitter.height)) +
      geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2, col=genes_info), arrow = arrow(length=unit(0.15,"cm"))) +
      geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") +
      geom_hline(yintercept = 0, lwd=0.2, linetype="dashed", col = "black") +
      geom_vline(xintercept = 0, lwd=0.2, linetype="dashed", col = "black") +
      ggtitle("PCA Loadings Plot") + scale_color_manual(values = c("#CC0000","orange")) +
      theme_bw() + my_theme + theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank()) 
      coord_equal(ratio = 1)
    plots <- list("all"=p,"genes"=p.genes)
    return(plots)  
  } else {
    return(p)
  }
  
}

# to be deprecated ----
plot_loadings <- function(pca_data, components, point_size = 0.5, assigned = NULL)
{
  if(is.list(pca_data)) {
    pca_load <- pca_data$pca$rotation
  } else {
    pca_load <- pca_data
  }
  
  toplot <- as.data.frame(pca_load[,components])
  if(!is.null(assigned)) {
    aidx <- match(assigned$variable, rownames(toplot))
    toplot$assigned <- assigned[aidx, 'princomp']
    nidx <- which(!toplot$assigned%in%components)
    toplot$assigned[nidx] <- "others"
    
    p0 <- ggplot(toplot, aes_string(x=components[1], y=components[2], col="assigned"))
  } else {
    p0 <- ggplot(toplot, aes_string(x=components[1], y=components[2]))
  }
  
  p <- p0 + geom_point(size=point_size, alpha=0.8) + 
    theme_bw() + my_theme + 
    geom_hline(yintercept = 0, lwd=0.2, linetype="dashed", col = "black") +
    geom_vline(xintercept = 0, lwd=0.2, linetype="dashed", col = "black") +
    ggtitle("PCA - Loading Plot") + scale_color_manual(values = c("#404040","#CC0000","orange")) +
    theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank()) 
  
  return(p)
}
assign_loadings <- function(pca_data, retained = NULL)
{
  if(is.list(pca_data)) {
    pca_load <- pca_data$pca$rotation
  } else {
    pca_load <- pca_data
  }
  
  if(!is.null(retained)) {
    if( is.numeric(retained) )  {
      if( length(retained)==1 )  retained <- 1:retained  
      retained <- paste0("PC",retained)
    }
    pca_load <- pca_load[,retained]
  }
  
  assigned <- apply(pca_load, 1, function(x) which.max(abs(x)))
  assigned <- paste0("PC",assigned)
  assigned <- data.frame("variable"    = rownames(pca_load)
                         , "princomp"  = assigned
                         , stringsAsFactors = F)
  return(assigned)
}
