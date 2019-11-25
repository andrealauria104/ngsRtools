# Count pattern from fasta sequence file ====
count_from_fasta <- function(fasta
                             , pattern = 'CG'
                             , outfile = NULL
                             , res = 'value' # value / fraction
                             # , ...
)
{
  require(Biostrings)
  require(parallel)
  
  message("[+] Count from fasta")
  message("  -- file: "   , fasta)
  message("  -- pattern: ", pattern)
  
  # sequence_list  <- sequence_to_list(fasta = fasta)
  #
  # pattern_counts <- mclapply(sequence_list
  #                            , vcountPattern
  #                            , pattern   = pattern
  #                            , mc.silent = T
  #                            , ... )
  fa <- readDNAStringSet(fasta)
  
  pattern_counts <- vcountPattern(pattern = pattern, subject = fa)
  names(pattern_counts) <- names(fa)
  
  pattern_counts <- reshape2::melt(pattern_counts)
  pattern_counts <- cbind.data.frame('sequence_id' = rownames(pattern_counts)
                                     ,  'pattern_counts' = pattern_counts[,'value'])
  if( res == 'fraction' ) {
    pattern_counts$pattern_counts <- pattern_counts$pattern_counts/width(fa)
  } else if ( res == 'value' ) {
    
  } else {
    stop(message('[!] Incorrect result type. (value/fraction)'))
  }
  
  if( !is.null(outfile) ) {
    if( grepl('\\.rds', outfile) ) {
      saveRDS(pattern_counts, file = outfile)
    } else {
      save(pattern_counts, file = outfile)
    }
  }
  
  return(pattern_counts)
}

plot_pattern <- function(pattern_counts, pattern
                         , ptitle = NULL)
{
  
  #source("theme_setting.R")
  
  pal <- pal_d3()(2)[2]
  
  if(any(grepl('GC content', pattern))) {
    lab <- pattern
  } else {
    lab <- paste0('Observed ',pattern, ' pattern')
  }
  
  bw <- function(x) (2 * IQR(x) / length(x)^(1/3)) # Freedman–Diaconis rule
  
  p <- ggplot(pattern_counts, aes(x=pattern_counts)) +
    geom_histogram(col='black', size = 0.25, alpha = 0.8, fill = pal, binwidth = bw) +
    theme_bw() + my_theme + xlab(lab)
  
  if(!is.null(ptitle)) {
    p <- p + ggtitle(ptitle) +
      theme(plot.title = element_text(size = 8, hjust = 0.5, face = 'bold'))
  }
  return(p)
}

# Match pattern from fasta sequence file ====
get_pattern_matching <- function(pattern_seq, fa, max.mismatch=1, ...)
{
  if(is.character(fa)) {
    message("[+] Reading from fasta: ", fa)
    fa <- readDNAStringSet(filepath = fa, format = "fasta")
  }
  pattern <- DNAStringSet(pattern_seq, use.names = T)
  
  matches <- vector(mode = 'list', length = length(pattern_seq))
  names(matches) <- names(pattern_seq)
  
  for(i in names(matches)) {
    message("[+] Matching sequence: ", pattern_seq[[i]])
    message(" -- max mismatch = ", max.mismatch)
    # match pattern
    match <- vmatchPattern(pattern   = pattern_seq[[i]]
                           , subject = fa
                           , with.indels  = F
                           , max.mismatch = max.mismatch
                           , ...
    )
    
    matches[[i]] <- as.data.frame(unlist(match))
    
    # get matched sequence
    idx      <- names(unlist(match))
    fa_match <- fa[match[idx]]
    seqs     <- unlist(lapply(idx, function(i) toString(fa_match[i])))
    names(seqs) <- idx
    
    ids <- unique(names(seqs)[which(duplicated(names(seqs)))])
    
    for(j in ids) {
      x <- unique(seqs[grep(j, names(seqs))])
      duplen <- nchar(x)/nchar(pattern_seq[[i]])
      if(duplen<2) {
        # print(duplen)
        print(j)
      }
      ss <- 1
      se <- nchar(pattern_seq[[i]])
      for(k in 1:duplen) {
        seqs[grep(j, names(seqs))][[k]] <- substr(x, ss, se)
        ss <- ss + nchar(pattern_seq[[i]])
        se <- se + nchar(pattern_seq[[i]])
      }
    }
    
    # retrieve chromosomal coordinates
    nm <- lapply(1:3, function(i) {
      y <- sapply(strsplit(idx, "[[:punct:]]"), "[[", i)
      if(i!=1) y <- as.numeric(y)
      return(y)}
    )
    
    # format output
    matches[[i]]$sequence <- seqs
    matches[[i]]$start <- nm[[2]] + matches[[i]]$start
    matches[[i]]$end   <- nm[[2]] + matches[[i]]$end
    matches[[i]] <- cbind.data.frame('chr' = nm[[1]], matches[[i]])
  }
  
  rm(match, fa)
  
  return(matches)
}
write_bed <- function(matches, bedfile)
{
  message("[+] Writing to: ", bedfile)
  write.table(matches[,c('chr','start','end','sequence')],
              quote = F      ,
              row.names = F  ,
              col.names = F  ,
              sep = "\t",
              file = bedfile)
}
wrapper_pattern_matching <- function(pattern_seq, fa, max.mismatch, outdir=NULL)
{
  matches        <- vector(mode = 'list', length = length(max.mismatch))
  names(matches) <- paste0('mismatch_', max.mismatch)
  for(m in max.mismatch) {
    matches[[m]] <- get_pattern_matching(pattern_seq = pattern_seq, fa = fa, max.mismatch = m)
    if(!is.null(outdir)) {
      # save bed --
      lapply(names(matches[[m]]), function(i) write_bed(matches[[m]][[i]], bedfile = paste0(outdir, "/",i,".mismatch_",m,".bed")))
    }
  }
  
  return(matches)
}
