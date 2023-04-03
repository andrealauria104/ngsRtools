plot_gene_expression <- function(y, genes
                                 , scale = F
                                 , lines_error_bar = T
                                 , gene_info_column = NULL
                                 , sample_name_column = NULL
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
    if(!any(colnames(y$samples)=="sample")) {
      if(!is.null(sample_name_column)) {
        y$samples$sample <- y$samples[,sample_name_column] 
      } else {
        y$samples$sample <- rownames(y$samples)
      }
    }
    if(!all(y$samples$sample==colnames(y))) {
      message("[!] Sample column names/metadata do not match.")
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
          m$sample <- apply(m[, summarize_by], 1, paste0, collapse = "_")
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
            pm <- ggplot(m, aes_string(x=summarize_by[1], y="av", group = summarize_by[2], col = summarize_by[2])) 
            if(lines_error_bar) pm <- pm + geom_errorbar(aes(ymax=av+se, ymin=av-se, group = summarize_by[2]), size=0.2, width=0.3,linetype="dashed", lwd = 0.25)
          }
          pm <- pm + geom_line(lwd=0.25) + geom_point() + 
            facet_wrap(as.formula(facet_formula), scales = "free_y", ncol = facet_ncol) + theme_bw() + my_theme_2 +
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
}

plotExpression <- function(y, gene
                           , experimental_info = NULL
                           , expression.unit = "CPM"
                           , group.by=NULL
                           , plot.type="bar"
                           , pal=NULL
                           , group_ord_idx=NULL
                           , ord_idx=NULL
                           , x_text_size=8
                           , show.points=T
                           , point.size=0.8
                           , point.width=0.1
                           , point.height=NULL
                           , col.width=0.8
                           , error.type=c("sd","se")
                           , show.names=F
                           , show.legend=T
                           , names.rot=0
                           , gene.text.face="bold") 
{
  plot.type <- match.arg(plot.type)
  error.type <- match.arg(error.type)
  
  if(class(y)=="DGEList") {
    expr <- y[[expression.unit]]
  } else {
    expr <- as.matrix(y)
  }
  
  if(is.null(experimental_info) && class(y)=="DGEList") {
    experimental_info <- y$samples
  }
  if(is.null(experimental_info)) {
    message("[!] Missing experimental_info.")
  } else {
    
    sample_col_idx <- grep("^sample$",colnames(experimental_info),ignore.case = T)
    if(length(sample_col_idx)!=0) {
      colnames(experimental_info)[sample_col_idx] <- "sample"  
    } else if(!is.null(rownames(experimental_info))) {
      experimental_info$sample <- rownames(experimental_info)
    }
    if(!identical(experimental_info$sample,colnames(expr))) {
      message("[!] Invalid experimental_info. Samples do not match.")
    } else {
      toplot <- reshape2::melt(expr[gene,,drop=F], varnames=c("gene","sample"))
      
      # browser()
      if(!is.null(group.by)) {
        message(" -- averaging over groups: ", group.by)
        colnames(experimental_info)[grep("^sample$",colnames(experimental_info),ignore.case = T)] <- "sample"
        toplot <- merge(toplot, experimental_info, by = "sample")
        
        toplot$group_by <- toplot[,group.by] # tmp variable to group by
        if(error.type=="sd") {
          toplot <- ddply(toplot, .(gene, group_by)
                          , mutate
                          , av = mean(value)
                          , sd = sd(value))
          
        } else if(error.type=="se") {
          toplot <- ddply(toplot, .(gene, group_by)
                          , mutate
                          , av = mean(value)
                          , sd = sd(value)/sqrt(length(value)))
        }
        toplot$group_by <- NULL # rm tmp variable
        
        if(!is.null(ord_idx)) toplot[,group.by] <- factor(toplot[,group.by], levels = ord_idx)
        if(is.null(pal)) pal <- ggsci::pal_d3()(length(unique(toplot[,group.by])))
        
        if(plot.type=="bar") {
          p <- ggplot(unique(toplot[,c(group.by,"av","sd","gene")]), aes_string(x=group.by,y="av", fill=group.by))+ 
            geom_col(lwd = 0.25, width=col.width, show.legend = show.legend) +
            geom_errorbar(aes(ymax=av+sd, ymin=av-sd), size=0.2, width=0.3,linetype="dashed", lwd = 0.25, position = position_dodge(width = 0.8)) +
            facet_wrap(~gene, scales = "free_y", ncol = 4) + theme_bw() + my_theme_2 +
            theme(legend.key.size = unit(4,'mm'),axis.title.x = element_blank(), strip.text = element_text(face = gene.text.face, size = 8)) +
            scale_fill_manual(values = pal) + ylab(paste0("average ",expression.unit))
          if(show.points) {
            p <- p + geom_jitter(data=toplot, aes_string(x=group.by,y="value"),width = point.width,height = point.height,size=point.size, show.legend=F) + 
              ylab(expression.unit)
          }
          if(show.names) {
            if(names.rot==45) {
              p <- p + theme(axis.text.x = element_text(size = x_text_size, angle=names.rot, hjust = 1, vjust = 1))
            } else if(names.rot==90){
              p <- p + theme(axis.text.x = element_text(size = x_text_size, angle=names.rot, hjust = 1, vjust = .5))
            }else {
              p <- p + theme(axis.text.x = element_text(size = x_text_size))
            }
          } else {
            p <- p + theme(axis.text.x = element_blank())
          }
        }
      } else {
        if(is.null(ord_idx)) {
          if(!is.null(group_ord_idx)) {
            ord_idx <- experimental_info$sample[order(experimental_info[,group_ord_idx])]
          } else {
            ord_idx <- colnames(expr)
          }
        }
        toplot$sample <- factor(toplot$sample, levels=ord_idx) 
        if(is.null(pal)) pal <- ggsci::pal_d3()(length(unique(toplot[,"sample"])))
        if(plot.type=="bar") {
          p <- ggplot(toplot, aes(x=sample, y=value, fill=sample))+ geom_col(lwd = 0.25,width=0.9,show.legend = F) +
            facet_wrap(~gene, scales = "free_y", ncol = 4) + theme_bw() + my_theme_2 +
            theme(legend.key.size = unit(4,'mm'),axis.title.x = element_blank(), axis.text.x = element_text(size=x_text_size, angle=45, hjust=1,vjust=1), strip.text = element_text(face = "plain", size = 8)) +
            scale_fill_manual(values = pal) + ylab(expression.unit)
        }
      }
      return(p)
    }
  }
}