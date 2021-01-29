plot_gene_expression <- function(y, genes
                                 , scale = F
                                 , lines_error_bar = T
                                 , gene_info_column = NULL
                                 , plot.type = 'lines'
                                 , expression.unit = "cpm"
                                 , summarize_by = NULL
                                 , group_by = 1
                                 , pal = NULL
                                 , precedence = NULL
                                 , facet_ncol = 4
                                 , facet_formula = "~gene_name"
                                 , show_text=T)  
{
  if(grepl("clean|correct|log",expression.unit,ignore.case = T)) {
    expr.unit.nm <- expression.unit
  } else {
    expr.unit.nm <- toupper(expression.unit)
  }

  if(class(y)!="DGEList") {
    message("[!] Invalid y, please provide DGEList object.")
  } else {
    if(is.data.frame(genes)) {
      if(!("gene_name"%in%colnames(genes))) {
        message("[!] Invalid genes data frame, it must contain \"gene_name\" column.")
      } else {
        gidx <- intersect(genes$gene_name, rownames(y[[expr.unit.nm]]))
      }
    } else if(is.character(genes)){
      gidx <- intersect(genes, rownames(y[[expr.unit.nm]]))
    } else if(is.list(genes) && !is.data.frame(genes)) {
      genes <- reshape2::melt(genes)
      colnames(genes) <- c("gene_name","gene_set")
      message(" -- reshaping gene set list, gene set column set as \"gene_set\"")
      gidx <- intersect(genes$gene_name, rownames(y[[expr.unit.nm]]))
    }
    if(length(gidx)!=0) {
      m <- y[[expr.unit.nm]][gidx,,drop=F]
      if(scale) {
        m <- t(scale(t(m)))
        expr.unit.nm <- paste0("scaled ", expr.unit.nm)
      }
      m <- reshape2::melt(m, varnames = c('gene_name','sample'))
      m <- merge(m, y$samples, by="sample")
      if(is.data.frame(genes)) m <- merge(m, genes, by="gene_name")
    } else {
      message("[!] Not detected gene list.")
    }
    
    if(!is.null(summarize_by)) {
      m <- ddply(m, c("gene_name",summarize_by)
                 , summarize
                 , av = mean(value)
                 , se = sd(value)/sqrt(length(value))
                 , .drop = F)
      if(is.data.frame(genes)) m <- merge(m, genes, by="gene_name")
      if (length(summarize_by) == 1) {
        m$sample <- m[,summarize_by]
      } else {
        m$sample <- apply(m$sample[, summarize_by], 1, paste0, collapse = "_")
      }
      if(!is.null(precedence)) m$sample <- factor(m$sample, levels = precedence)
      if(plot.type=="bars") {
        if(is.null(pal)) pal <- ggsci::pal_d3()(length(unique(m$sample)))
        if(show_text) {
          axis_text_x <- element_text(angle=45,hjust=1,vjust=1,size=6)
        } else {
          axis_text_x <- element_blank()
        }
        pm <- ggplot(m, aes(x=sample, y=av, fill=sample)) + geom_col(lwd = 0.25) +
          geom_errorbar(aes(ymax=av+se, ymin=av-se), size=0.2, width=0.3,linetype="dashed", lwd = 0.25, position = position_dodge(width = 0.8)) +
          facet_wrap(as.formula(facet_formula), scales = "free_y", ncol = facet_ncol) + theme_bw() + my_theme_2 +
          theme(axis.title.x = element_blank(), axis.text.x = axis_text_x, strip.text = element_text(face = "plain", size = 8)) +
          scale_fill_manual(values = pal) + ylab(paste0("average ",expr.unit.nm))
      } else if(plot.type=="lines") {
        if(length(summarize_by) == 1) {
          if(group_by=="gene_name") {
            col_by <- "gene_name"
          } else {
            col_by <- "sample"
          }
          if(is.null(pal)) pal <- ggsci::pal_d3()(length(unique(m[,col_by])))
          pm <- ggplot(m, aes_string(x="sample", y="av", col = col_by,group = group_by))
          if(lines_error_bar) pm <- pm + geom_errorbar(aes(ymax=av+se, ymin=av-se, group = summarize_by[2]), size=0.2, width=0.3,linetype="dashed", lwd = 0.25)
        } else {
          if(is.null(pal)) pal <- ggsci::pal_d3()(length(unique(m[,summarize_by[2]])))
          pm <- ggplot(m, aes_string(x=summarize_by[1], y=av, group = summarize_by[2], col = summarize_by[2])) 
          if(lines_error_bar) pm <- pm + geom_errorbar(aes(ymax=av+se, ymin=av-se, group = summarize_by[2]), size=0.2, width=0.3,linetype="dashed", lwd = 0.25)
        }
        pm <- pm + geom_line(lwd=0.25) + geom_point() + 
          facet_wrap(as.formula(facet_formula), scales = "free_y", ncol = 6) + theme_bw() + my_theme_2 +
          scale_color_manual(values = pal) + ylab(paste0("average ",expr.unit.nm))
      }} else {
        if(plot.type=="bars") {
          if(is.null(pal)) pal <- ggsci::pal_d3()(length(unique(m$sample)))
          pm <- ggplot(m, aes(x=sample, y=value, fill=sample)) + geom_col(lwd = 0.25) +
            facet_wrap(as.formula(facet_formula), scales = "free_y", ncol = facet_ncol) + theme_bw() + my_theme_2 +
            theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=6), strip.text = element_text(face = "plain", size = 8)) +
            scale_fill_manual(values = pal) + ylab(expr.unit.nm)
        } else if(plot.type=="lines") {
          message("[!] Invalid plot type, provide \"summarize_by\" for line plotting.")
        }}
    return(pm)
  }
}