# Differential Expression Analysis ====
#' Differential Expression Analysis with edgeR
#'
#' Implementation of statistical models and tests 
#' for differential expresion analysis using edgeR
#' software package
#' @param y count matrix (or data.frame) or object of class DGEList
#' @param experimental_info data.frame containing informations about the experimental conditions for each sample
#' @param group conditions to test
#' @param method statistcal test to be performed (exact, lrt, qlf)
#' @return results of the test (data.frame)
#'
#' @examples
#' \dontrun{
#'     dea <- calculateDiffExprEdgeR(y, experimental_info
#'     , group = experimental_info$group)
#' }
#' 
#'
#' @export
calculateDiffExprEdgeR <- function(y
                                   , experimental_info = NULL
                                   , gene_info = NULL
                                   , group  = NULL
                                   , reference = NULL # Control group
                                   , filter = T
                                   , filter.expr.th = 1
                                   , filter.sample.th = 2
                                   , expression.unit = c("cpm","rpkm")
                                   , tlen = NULL
                                   , method = c("qlf","lrt","exact")
                                   , robust.dispersion = TRUE
                                   , robust.ql = FALSE
                                   , anovalike = F # Activate ANOVA-like test for any difference (alternative to cf,contrast)
                                   , design = NULL # Custom design matrix
                                   , cf = NULL # Testing coefficient (default: last column in design matrix)
                                   , contrast = NULL # Contrast matrix - character/numeric vector
                                   , return.y  = F # Return y? For backward compatibility, better set return.res.type
                                   , return.res.type = c("de","y","all") # What to return
                                   ) {
  
  
  # Compute differential expression with edgeR
  expression.unit <- match.arg(expression.unit)
  method <- match.arg(method)
  return.res.type <- match.arg(return.res.type)
  
  # Internals ---
  get_contrast <- function(contrast, design)
  {
    if(is.character(contrast) & length(contrast) == 3) {
      cnm <- paste0(contrast[1], contrast[-1], collapse = "-")
      message(" -- Contrast: ", cnm)
      mcontrast <- limma::makeContrasts(contrasts = cnm, levels = design)
      print(mcontrast)
      return(mcontrast)  
    } else {
      print(contrast)
      return(contrast)
    }
    
  }
  print_coefficients <- function(cf) {
    if(is.character(cf)) {
      nmcf <- cf
    } else {
      nmcf <- colnames(design)[cf]
    }
    if(length(nmcf)>1) nmcf <- paste0(nmcf,collapse = "-")
    message(" -- Testing coefficient: ",nmcf)
  }
  
  message("[*] Run edgeR for Differential Expression Analysis")
  
  if(is.matrix(y) | is.data.frame(y)) {
    y <- processRNAseqEdgeR(m = y
                            , experimental_info = experimental_info
                            , gene_info = gene_info
                            , group = group
                            , reference = reference
                            , filter = filter
                            , filter.expr.th = filter.expr.th
                            , filter.sample.th = filter.sample.th
                            , expression.unit = expression.unit
                            , tlen = tlen
                            , return.normalized.expr = F)
  }
  
  if(is.null(design)) {
    # Standard design matrix 
    # First column is control ( = reference in group)
    message(" -- Standard design matrix: ~group")
    if(!any(grepl('group', colnames(y$samples)))) {
      if(!is.null(group)) {
        if(length(group)==1) {
          y$samples$group <- y$samples[,group]  
        } else {
          y$samples$group <- apply(y$samples[,group]  , 1, paste0, collapse = "_")
        }
        y$samples$group <- factor(y$samples$group)  
        if(!is.null(reference)) y$samples$group <- relevel(y$samples$group, ref = reference) 
      } else {
        stop(message("[-] group not available in metadata. Please, specify a valid design matrix."))
      }
    }
    design <- model.matrix(~group, data = y$samples)
  } else {
    # Customized design matrix
    # Passed as string formula or model.matrix
    message(" -- Parsing design matrix")
    if(is.character(design)) {
      message("    * Formula: ",design)
      design <- model.matrix(as.formula(design), data = y$samples)
    }
  }
  
  if(!identical(rownames(design), colnames(y))) {
    message("[!] Setting rownames in design matrix equal to sample names.")
    rownames(design) <- colnames(y)
  }
  
  if(is.null(cf) & is.null(contrast)) {
    if(anovalike) {
      # Test all columns in design matrix (ANOVA-like)
      # parametrize model wrt control ( = reference in group)
      cf  <- 2:ncol(design)
    } else {
      # Default: last column in design matrix
      cf  <- ncol(design)
    }
  }
  
  message(" -- Estimating dispersion, robust = ", robust.dispersion)
  y <- edgeR::estimateDisp(y, design, robust=robust.dispersion)
  
  message(" -- Testing differential expression, method: ", method)
  if(length(cf)>1) {
    message(" -- ANOVA-like for multiple group comparison") 
  }
  message(" -- Design matrix: ")
  print(design)
  
  if(method=="exact") {
    # Exact test
    et <- edgeR::exactTest(y)
    de <- edgeR::topTags(et, n = Inf)
    
  } else if(method=="lrt") {
    # Likelihood-ratio test
    fit <- edgeR::glmFit(y, design)
    if(is.null(contrast)){
      # Standard comparison
      print_coefficients(cf)
      lrt <- edgeR::glmLRT(fit, coef=cf)
    } else {
      # GLM with contrasts
      mcontrast <- get_contrast(contrast = contrast, design = design)
      lrt <- edgeR::glmLRT(fit, contrast = mcontrast)
    }
    de <- edgeR::topTags(lrt, n = Inf)
    
  } else if(method=="qlf") {
    # Quasi-likelihood F-test
    message(" -- glmQLFit, robust QL = ", robust.ql)
    fit <- edgeR::glmQLFit(y, design, robust=robust.ql)
    if(is.null(contrast)){
      # Standard comparison
      print_coefficients(cf)
      qlf <- edgeR::glmQLFTest(fit, coef=cf)
    } else {
      # GLM with contrasts
      mcontrast <- get_contrast(contrast = contrast, design = design)
      qlf <- edgeR::glmQLFTest(fit, contrast = mcontrast)
    }
    de <- edgeR::topTags(qlf, n = Inf)
    
  } else {
    message("[-] Test method not available. Please, provide a valid one ( exact / lrt / qlf ).")
  }
  
  if(return.y || return.res.type=="y") {
    
    res <- list("y" = y, "de" = de)
    return(res)
    
  } else if(return.res.type=="all" &&  method=="qlf") {
    
    res <- list("y"=y, "de"=de, "fit"=fit, "qlf"=qlf)
    return(res)
    
  } else if(return.res.type=="all" &&  method=="lrt") {
    
    res <- list("y"=y, "de"=de, "fit"=fit, "lrt"=lrt)
    return(res)
    
  } else if(return.res.type=="all" &&  method=="exact") {
    
    res <- list("y"=y, "de"=de, "et"=et)
    return(res)
    
  } else if(return.res.type=="de") {
    
    return(de)
    
  }
}

saveXLSresEdgeR <- function(res, outfile, name, force=T) {
  # set java memory
  options(java.parameters = "-Xmx8000m")
  # java garbage collector
  # jgc <- function()
  # {
  #   .jcall("java/lang/System", method = "gc")
  # }
  
  message("[+] Saving results to file: ", outfile)
  
  if( is.list(res) & !is.data.frame(res)) {
    for( i in 1:length(res) ) {
      
      gc(verbose = F)
      # jgc()
      tox <- res[[i]]
      if(i==1 & file.exists(outfile) & force==F) {
        ap <- T
        warning(" -- appending to existing file.")
      } else {
        ap <- i>1
      }
      xlsx::write.xlsx2(tox, file = outfile, sheetName = names(res)[i], append = ap, row.names = F)
    }
  } else {
    # ap  <- F
    if(!is.data.frame(res)) res <- as.data.frame(res)
    tox <- res
    if(missing(name)) name <- "edgeR"
    if(file.exists(outfile) & force==F) {
      ap <- T
      warning(" -- appending to existing file.")
    } else {
      ap <- F
    }
    xlsx::write.xlsx2(tox, file = outfile, sheetName = name, append = ap, row.names = F)
  }
  
}

saveXLSresEdgeR2 <- function(res, outfile
                             , name = "edgeR"
                             , force = TRUE
                             , keepNA = TRUE)
{
  message("[+] Saving results to file: ", outfile)
  if(is.list(res) & !is.data.frame(res)) {
    
    if(file.exists(outfile) & force==F) {
      warning(" -- appending to existing file.")
      wb <- openxlsx::loadWorkbook(outfile)
    } else {
      wb <- openxlsx::createWorkbook()
    }
    
    for(i in names(res)) {
      openxlsx::addWorksheet(wb = wb, sheetName = i)
      openxlsx::writeData(wb = wb
                          , sheet = i
                          , x = res[[i]]
                          , keepNA = keepNA
                          , colNames = TRUE
                          , rowNames = FALSE)
    }
    
    openxlsx::saveWorkbook(wb = wb
                           , file = outfile
                           , overwrite = TRUE)
    
  } else {
    
    if(file.exists(outfile) & force==F) {
      overwrite <- T
      warning(" -- appending to existing file.")
    } else {
      overwrite <- F
    }
    openxlsx::write.xlsx(x = res
                         , file = outfile
                         , overwrite = TRUE
                         , keepNA = keepNA
                         , sheetName = name
                         , colNames = TRUE
                         , rowNames = FALSE)
    
  }
}

getDEgenes <- function(x, fdrTh=0.1, fcTh=0.5, lcpmTh=NULL) {
  # Get DE genes satisfying thresholds
  message("[+] Get differentially expressed genes, thresholds:")
  message(" - logFC = ", fcTh)
  message(" - FDR = ", fdrTh)
  message(" - logCPM = ", lcpmTh)
  
  if(!is.data.frame(x)) x <- x$table
  
  if( !is.null(fdrTh) ) {
    x <- x[ x[,'FDR'] <= fdrTh , ]
  }
  
  if( !is.null(lcpmTh) ) {
    x <- x[ x[,'logCPM'] >= lcpmTh , ]
  }
  
  if(nrow(x) == 0) {
    stop(message("[-] No differentially expressed genes with provided FDR"))
  }
  
  if( !is.null(fcTh) ) {
    fcidx <- grep("logFC", colnames(x))
    
    if( length(fcidx)==1 ) {
      x <- x[ abs(x[,fcidx]) >= fcTh , ]
    } else {
      x <- x[ rowSums( abs(x[,fcidx]) >= fcTh) >= 1,]
    }
    
    if(nrow(x) == 0) {
      stop(message("[-] No differentially exxpressed genes with provided logFC"))
    }
  }
  
  return(x)
  
}

getDEgsigned <- function(deg, signed=(-1))
{
  fcidx <- grep("logFC", colnames(deg))
  x     <- deg[,fcidx]
  
  max.idx <- apply(x, 1, function(m) which.max(abs(m)))
  
  y <- vector(mode = 'numeric', length = nrow(x))
  names(y) <- names(max.idx)
  
  for(i in names(max.idx)) {
    y[i] <- x[i, max.idx[i]]
  }
  
  gene.idx <- names(y[which(sign(y)==signed)])
  
  return(x[gene.idx,])
}

plotDiffExprRes <- function(de
                            , type = c("volcano","MA")
                            , lfcTh = 1
                            , pvTh = 0.05
                            , top = 5
                            , gtitle = NULL
                            , pal = NULL
                            , label_selection = NULL
                            , corrected_pval = T
                            , point.size = 0.5)
{
  type <- match.arg(type)
  
  # plot differential expression results for edgeR, DESeq2, limma
  if(is.null(pal)) pal <- c('#004C99','#404040','#CC0000')
  
  fcIdx <- grep("logFC|log2FoldChange", colnames(de), value = T)
  
  if(corrected_pval) {
    pvIdx <- grep("adj.P.Val|FDR|padj", colnames(de), value = T)
    volcano_lab <- "adjusted P-value"
  } else {
    pvIdx <- grep("P.Value|PValue|pvalue", colnames(de), value = T)
    volcano_lab <- "P-value"
  }
  message(" -- significance threshold on: ", volcano_lab)
  
  cIdx <- grep("baseMean|logCPM|AveExpr", colnames(de), value = T)
  
  if(length(fcIdx)!=1 || length(pvIdx)!=1 || length(cIdx)!=1) {
    message("[!] Invalid format:")
    if(length(fcIdx)!=1) {
      message(" -- more than 1 valid Fold Change column")
      print(fcIdx)
    }
    if(length(pvIdx)!=1) {
      message(" -- more than 1 valid PValue column")
      print(pvIdx) 
    }
    if(length(cIdx)!=1) {
      message(" -- more than 1 valid average expression column")
      print(cIdx)
    }
  } else {
    tmp   <- de[,c(cIdx,fcIdx,pvIdx)]
    colnames(tmp) <- c("meanc","lfc","padj")
    tmp$status <- "none"
    tmp$status[which(tmp$lfc>=lfcTh & tmp$padj<=pvTh)] <- "Up-regulated"
    tmp$status[which(tmp$lfc<=(-lfcTh) & tmp$padj<=pvTh)] <- "Down-regulated"
    tmp$status <- factor(tmp$status)
    
    topUP <- tmp[tmp[,'status']=="Up-regulated",]
    nup <- nrow(topUP)
    topDW <- tmp[tmp[,'status']=="Down-regulated",]
    ndw <- nrow(topDW)
    
    if(!is.null(label_selection)) {
      tmp$lab <- rownames(tmp)
      tmp$lab[-which(rownames(tmp)%in%label_selection)] <- ""
      
    } else if(!is.null(top)) {
      topUP <- rownames(topUP[1:top,])
      topDW <- rownames(topDW[1:top,])
      
      tmp$lab <- rownames(tmp)
      tmp$lab[-which(rownames(tmp)%in%c(topDW, topUP))] <- ""
    } else {
      tmp$lab <- ""
    }
    
    if(type=="volcano") {
      
      message(" -- plot volcano, thresholds:")
      message(" - |logFC| = ", lfcTh)
      message(" - pvalue = ", pvTh)
      
      p <- ggplot(tmp, aes(x=lfc, y=-log10(padj),col=status, label=lab)) + geom_point(size=point.size, alpha=0.6) + theme_bw() +
        geom_hline(yintercept = -log10(pvTh), linetype = 'dashed', lwd = 0.25) + 
        geom_vline(xintercept = lfcTh, linetype = 'dashed', lwd = 0.25) +
        geom_vline(xintercept = -lfcTh, linetype = 'dashed', lwd = 0.25) +
        xlab("logFC") + ylab(volcano_lab) +
        theme_bw() + my_theme_2 + ggtitle(gtitle) +
        scale_color_manual(values = pal) +
        ggrepel::geom_text_repel(show.legend = F, size = 2, segment.size = 0.1) +
        geom_label(
          data    = subset(tmp, status=="Up-regulated"),
          mapping = aes(x = max(tmp$lfc)-.5, y = 1, label = nup, col = status),
          size = 2,
          show.legend = F
        ) +
        geom_label(
          data    = subset(tmp, status=="Down-regulated"),
          mapping = aes(x = min(tmp$lfc)+.5, y = 1, label = ndw, col = status),
          size = 2,
          show.legend = F
        ) +
        guides(col = guide_legend(nrow = 2))
    } else if(type=="MA") {
      
      message(" -- plot MA, thresholds:")
      message(" - |logFC| = ", lfcTh)
      message(" - pvalue = ", pvTh)
      
      if(cIdx=="baseMean") tmp$meanc <- log2(tmp$meanc)
      p0 <- ggplot(data = tmp, aes(x=meanc, y=lfc, col=status, label=lab)) + geom_point(data=subset(tmp, status=="none"), aes(x=meanc, y=lfc, col=status, label=lab), size=point.size, alpha = 0.5)
      p1 <- geom_point(data= subset(tmp, status!="none"), aes(x=meanc, y=lfc, col=status, label=lab), size=point.size, alpha = 0.9)
      
      p <- p0 + p1 + theme_bw() +
        geom_hline(yintercept = 0, linetype = 'dashed', lwd = 0.25) + 
        ggrepel::geom_text_repel(show.legend = F, size = 2, segment.size = 0.1) +
        ylab("logFC") + xlab("mean normalized counts") +
        theme_bw() + my_theme_2 + ggtitle(gtitle) +
        scale_color_manual(values = pal) +
        guides(col = guide_legend(nrow = 2))
    }
    return(p)
  }

}

build_model_matrix <- function(group, ref, formula_string)
{
  # build model --
  if(is.character(group) & !missing(ref)) {
    # parametrize with "ref" as baseline 
    group <- factor(group)
    group <- relevel(group, ref = ref)
    mod   <- model.matrix(~group)
    colnames(mod)[-1] <- paste0(levels(group)[-1], "vs", levels(group)[1])
  } else {
    if(missing(formula_string)) {
      stop(message("[!] Please, provide parametrization formula (as charachter)"))
    }
    fm  <- as.formula(formula_string)
    mod <- model.matrix(fm, data = group)
  }
  
  return(mod)
}

prepare_pairwise_contrasts <- function(design, group, exclude = NULL)
{
  if(missing(design)) {
    group  <- as.factor(group)  
    design <- model.matrix(~group)
  }
  
  pairwise  <- apply(combn(rev(colnames(design)), 2),2, paste, collapse = "-")
  pairwise  <- gsub("\\-\\(Intercept\\)","", pairwise)
  if(!is.null(exclude)) {
    pairwise <- pairwise[grep(exclude, pairwise, invert = T)]
  }
  contrasts <- limma::makeContrasts(contrasts = pairwise, levels = design)
  
  return(list("contrasts" = contrasts, "design" = design))
}

fix_de_contrast_names <- function(de, ref)
{
  idx_nm <- grep("-", names(de))
  tmp_nm <- names(de)[idx_nm]
  names(de)[idx_nm] <- gsub("-","vs",gsub(paste0("vs",ref),"",tmp_nm))
  
  return(de)
}

analyze_de_contrasts_v2 <- function(y, test_contrasts
                                 , fdrTh   = 0.05
                                 , fcTh    = 1
                                 , lcpmTh  = 0.5
                                 , fixnm_ref = NULL
                                 , ...)
{
  
  de <- lapply(test_contrasts, function(cnt) {
    
    message("Contrast: ", paste0(cnt[1],cnt[-1], collapse = '-'))
    
    x <- calculateDiffExprEdgeR(y = y
                                , contrast = cnt
                                , ... )
    
    sig <- tryCatch(expr = getDEgenes(x, fdrTh = fdrTh, fcTh = fcTh, lcpmTh = lcpmTh)
                  , error = function(e) {message(e);return(NA)})
    res <- list("table" = x$table, "sig" = sig)
    return(res)
  })
  
  if(!is.null(fixnm_ref)) de <- fix_de_contrast_names(de, ref = fixnm_ref)
  
  return(de)
}

create_fc_matrix_contrasts <- function(de_contrasts)
{
  gidx <- lapply(de_contrasts, function(i) i$sig$genes)
  gidx <- unique(unlist(gidx))
  
  fc <- lapply(de_contrasts, function(i) i$table[gidx,"logFC"])
  fc <- do.call(cbind, fc)
  rownames(fc) <- gidx
  
  return(fc)
}

compute_fold_change <- function(x, y, pseudocount, logscale = F) 
{
  if(logscale) {
    fc <- x-y
  } else {
    fc <- log2((x+pseudocount)/(y+pseudocount))  
  }
  return(fc)
} 

analyze_relative_expression <- function(genes, dea, time_0 = "T0h", facet_formula = NULL)
{
  get_expression_fc <- function(dea, genes)
  {
    summarize_fc <- function(x, direction)
    {
      nx <- nrow(x)
      x <- as.data.frame(t(apply(x, 2, data_summary)))
      x$time <- rownames(x)
      x$direction <- direction
      x$n <- nx
      
      return(x)
    }
    
    if(!missing(genes)) {
      sidx    <- intersect(genes, dea$sig$genes)
      dea$sig <- dea$sig[sidx,]
      dea$table <- dea$table[genes,]
    }
    up <- getDEgsigned(dea$sig, signed = 1)
    gup <- rownames(up)
    up <- summarize_fc(x = up, direction = "up")
    dw <- getDEgsigned(dea$sig, signed = (-1))
    gdw <- rownames(dw)
    dw <- summarize_fc(x = dw, direction = "dw")
    
    cgenes <- setdiff(dea$table$genes, dea$sig$genes)
    const  <- dea$table[cgenes, grep("FC", colnames(dea$table))]
    gconst <- rownames(const)
    const <- summarize_fc(x = const, direction = "const")
    
    toplot <- rbind.data.frame(up, dw, const)  
    toplot <- plyr::unrowname(toplot)
    toplot$time <- gsub("logFC.","", toplot$time)
    toplot$time <- sapply(strsplit(toplot$time, "_"),"[[",1)
    
    t0 <- data.frame("y" = rep(0,3)
                     , "ymin" = rep(0,3)
                     , "ymax" = rep(0,3)
                     , "time" = rep(time_0,3)
                     , "direction" = c("up","dw","const")
                     , "n" = c(unique(up$n),unique(dw$n),unique(const$n)))
    
    toplot <- rbind.data.frame(toplot, t0)  
    
    return(list("toplot"=toplot, "gup" = gup, "gdw" = gdw, "gconst" = gconst))
  }
  wrap_get_efc <- function(dea, genes)
  {
    efc <- lapply(dea, function(i) 
    {
      efc <- get_expression_fc(i, genes)
      
    })
    
    reg <- lapply(efc, function(i) {
      i$toplot <- NULL
      i
    } )
    
    efc <- lapply(efc, "[[", "toplot")
    efc <- do.call(rbind, efc)
    efc$cond <- gsub("\\.+\\d+", "",rownames(efc))
    
    return(list("efc"=efc,"reg"=reg))
  }
  plot_expression_fc <- function(efc, facet_formula = NULL)
  {
    if(is.list(efc) & !is.data.frame(efc)) {
      efc <- do.call(rbind, efc) 
      efc$group <- sapply(strsplit(rownames(efc), "[.]"),"[[",1)
      if(is.null(facet_formula)) {
        p0 <- ggplot(efc, aes(x=time, y=y)) + facet_grid(group~cond) 
      } else {
        p0 <- ggplot(efc, aes(x=time, y=y)) + facet_wrap(as.formula(facet_formula), ncol = 3) 
      }
    } else {
      p0 <- ggplot(efc, aes(x=time, y=y)) + facet_grid(~cond) 
    }
    
    p <- p0 +
      theme_bw() + my_theme + ylab("average logFC") +
      geom_line(aes(group = direction, col=direction), size = 0.5) +
      geom_label(
        data    = subset(efc, direction=="up"),
        mapping = aes(x = 1, y = 3, label = n, col = direction),
        size = 3
      ) +
      geom_label(
        data    = subset(efc, direction=="dw"),
        mapping = aes(x = 1, y = -2, label = n, col = direction),
        size = 3
      ) +
      geom_ribbon(aes(ymin=ymin, ymax=ymax, group=direction, fill=direction), alpha = 0.2) +
      scale_fill_manual(values = c("lightgrey","blue","red")) +
      scale_color_manual(values = c("lightgrey","blue","red")) +
      theme(panel.grid = element_blank()
            , axis.title.x = element_blank()
            , axis.text.x = element_text(angle = 90, hjust = 1)
            , strip.background = element_blank()
            , legend.position =  "none"
            # , strip.text = element_text(face="bold", size = 8)
      ) 
    
    return(p)
  }
  
  if(is.list(genes)) {
    tmp <- lapply(genes, wrap_get_efc, dea = dea)
    efc <- lapply(tmp, "[[", "efc")
    reg <- lapply(tmp, "[[", "reg")
    rm(tmp)
    efc <- list("efc"=efc,"reg"=reg)
  } else {
    efc <- wrap_get_efc(dea = dea, genes = genes)
  }
  
  p <- plot_expression_fc(efc$efc, facet_formula = facet_formula)
  
  return(list("p" = p, "data" = efc))
} 

# DESeq2 ---
calculateDiffExprDESeq2 <- function(counts
                                    , info_analysis
                                    , design_formula 
                                    , contrast_tests = NULL
                                    , fcth           = NULL
                                    , pvth           = NULL
                                    , outfile        = NULL
                                    , ...) 
{
  message("[*] Run DESeq2 for Differential Expression Analysis")
  m           <- counts[,rownames(info_analysis)]
  m           <- apply(m, 2, as.integer)
  rownames(m) <- rownames(counts)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m,
                                        colData   = info_analysis,
                                        design    = as.formula(design_formula))
  dds <- DESeq2::DESeq(dds, ...)
  
  if(!is.null(contrast_tests)) {
    res <- list()
    for(i in contrast_tests) {
      k <- paste0(i[1],"_",i[2],"_vs_",i[3])
      if(!is.null(fcth) & !is.null(pvth)) {
        res[[k]] <- list()
        res[[k]]$result <- DESeq2::results(dds, contrast = i)
        res[[k]]$sig <- subset(DESeq2::results(dds, contrast = i), abs(log2FoldChange)>=fcth & padj <= pvth)
        res[[k]]$sig <- cbind.data.frame("genes"=rownames(res[[k]]$sig),res[[k]]$sig)
        res[[k]]$sig <- res[[k]]$sig[order(res[[k]]$sig$padj, decreasing = F),]
      } else {
        res[[k]] <- DESeq2::results(dds, contrast = i)
      }
    }
    if(!is.null(outfile)) {
      saveXLSresEdgeR(lapply(res,"[[","sig"), outfile = outfile, force = T)
    }
    
    return(list("dds"=dds,"res"=res))
  } else {
    return(dds)
  }
}
