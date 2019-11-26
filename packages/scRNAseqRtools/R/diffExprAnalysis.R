# Differential Expression ----
plot_cmarker_expression <- function(ebs, genes, pal, structure = "Cluster",log = T, assay.type = NULL, point.size = 0.5, plot.type = 'boxplot', filter = F)
{
  if(is.null(assay.type)) {
    if(log) {
      marker <- logcounts(ebs)[genes,]
    } else {
      marker <- normcounts(ebs)[genes,]
    }
  } else {
    marker <- assay(ebs, assay.type)[genes,]
  }
  
  marker <- reshape2::melt(marker)
  colnames(marker) <- c("name","cell","expression")
  marker$cluster <- colData(ebs)[marker$cell,grep(structure, colnames(colData(ebs)))]
  if(filter) marker <- marker[marker$expression>0,]
  if(plot.type=='boxplot') {
    dodge <- position_dodge(width = 0.8)
    p0 <- ggplot(marker, aes(x=cluster, y=expression, col = cluster)) + geom_boxplot(lwd = 0.25, outlier.size = 0.1) + geom_jitter(size = point.size)
  } else if(plot.type=='violin') {
    dodge <- position_dodge(width = 0.8)
    p0 <- ggplot(marker, aes(x=cluster, y=expression, col = cluster)) +
      geom_violin(trim = T, position = dodge) + geom_boxplot(width=0.08, outlier.color = NA, position = dodge) +
      stat_summary(fun.y=median, geom="point", size=0.5, color="black", position = dodge)
  }
  p <- p0 +
    facet_wrap(~name, ncol = 3) +
    theme_bw() + my_theme + scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) + xlab(structure) + ylab("Expression")  +
    theme(panel.grid = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1
          , strip.background = element_blank()
          , strip.text = element_text(size = 8, face = "bold"))
  
  return(p)
}

run_differrential_expression <- function(i, cds, structure, ntop, direction = "up", rank.var = "pval", ...)
{
  pData(cds)$category <- pData(cds)[,structure] == i
  cgroup_i <- rownames(pData(cds))[which(pData(cds)$category==T)]
  cgroup_rest <- rownames(pData(cds))[which(pData(cds)$category==F)]
  mean_1 <- rowMeans(exprs(cds)[,cgroup_i])
  mean_2 <- rowMeans(exprs(cds)[,cgroup_rest])
  
  test  <- differentialGeneTest(cds, fullModelFormulaStr = "~ category", ...)
  test$mean_cgroup_i <- mean_1[rownames(test)]
  test$mean_cgroup_rest <- mean_2[rownames(test)]
  test$logFC <- log2(test$mean_cgroup_i+1) - log2(test$mean_cgroup_rest+1)
  
  if(direction=="up") {
    results <- subset(test, mean_cgroup_i>mean_cgroup_rest)
  } else if(direction=="dw") {
    results <- subset(test, mean_cgroup_i<mean_cgroup_rest)
  } else if(direction=="all") {
    results <- test
  }
  if(!is.null(ntop)) results <- head(results[order(results[,rank.var], decreasing = rank.var=="logFC"),], n = ntop)
  
  return(results)
}

find_marker_genes <- function(cds, structure = "Cluster", ntop = 50, names_only = T, ...)
{
  structures <- levels(pData(cds)[,structure])
  structure_markers <- lapply(structures, run_differrential_expression
                              , cds  = cds
                              , ntop = ntop
                              , structure = structure
                              , ...)
  names(structure_markers) <- structures
  
  if((!is.null(ntop) || ntop > 0) && names_only == T) {
    top <- lapply(structure_markers, "[[", "feature_symbol")
    top <- as.data.frame(do.call(cbind, top))
    
    return(top)
  } else {
    return(structure_markers)
  }
  
  
}