# Ranking metric ====
getRankingMetric <- function(x, dea.tool="edgeR", method = "logFC_PValue")
{
  if(dea.tool=="edgeR") {
  	if(method=="logFC_PValue") {
  		message(" -- ranking metric: logFC * (-10*log(PValue)) [default]")
  		rank <- with(x, logFC * (-10*log(PValue)))
  	} else if(method=="logFC") {
  		message(" -- ranking metric: logFC")
  		warning("[!] By ranking with logFC only, you are not using the statistical information resulting from DE analysis.")
  		rank <- x[,"logFC"]
  	} else if(method=="PValue") {
  		message(" -- ranking metric: sign(logFC) * (-10*log(PValue))")
  		rank <- with(x, sign(logFC) * (-10*log(PValue)))
  	} else if(method=="FDR") {
  		message(" -- ranking metric: sign(logFC) * (-10*log(FDR))")
  		rank <- with(x, sign(logFC) * (-10*log(FDR)))
  	}else if(method=="lrt") {
  		message(" -- ranking metric: sign(logFC) * LR")
  		rank <- with(x, sign(logFC) * LR)  
  	} else if(method=="qlf") {
  		message(" -- ranking metric: sign(logFC) * F")
  		rank <- sign(x[,"logFC"]) * x[,"F"]
  	} else {
  		message("[!] Invalid method for edgeR; provide one of: logFC_PValue, logFC, PValue, FDR, lrt, qlf")
  	}
  }

  if(exists("rank",inherits = F)) {
    names(rank) <- x[,1]
    rank <- rank[order(rank,decreasing = F)]
    return(rank)
  }
  
}

# Fast GSEA ====
plot_fgseaRes <- function(fgseaRes,topPathways=NULL) 
{
  toplot <- as.data.frame(fgseaRes[,c(1,5)])
  if(!is.null(topPathways))  toplot <- subset(toplot, pathway%in%topPathways)
  toplot$direction <- factor(sign(toplot$NES))
  toplot$pathway <- gsub("_"," ",toplot$pathway)
  toplot$pathway <- gsub("MM ","",toplot$pathway)
  toplot$pathway <- gsub("(KEGG|REACTOME)","\\1:",toplot$pathway)
  toplot$pathway <- factor(toplot$pathway, levels=toplot$pathway[order(toplot$NES)])
  ggplot(toplot, aes(x=pathway, y=NES, fill=direction)) + geom_col() + coord_flip() +
    theme_bw() + my_theme + scale_fill_manual(values=c("#004C99", "#CC0000"))
  
  
}
get_top_fgseaRes <- function(fgseaRes, ntop=10) 
{
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=ntop), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=ntop), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
}