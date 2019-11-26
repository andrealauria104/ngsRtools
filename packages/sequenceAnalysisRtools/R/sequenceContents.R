# GC content ====
compute_gc_content <- function(sequence)
{
  n <- sapply(Biostrings::DNA_BASES, Biostrings::vcountPattern, subject=sequence)
  rownames(n) <- names(sequence)
  gc_content  <- round((n[,'G']+n[,'C'])/rowSums(n),4)
  
  return(gc_content)
}

get_gc_from_fasta <- function(fasta
                              , outfile = NULL
                              # , ...
)
{
  message("[+] Compute GC content from fasta")
  message("  -- file: "   , fasta)
  
  sequence   <- Biostrings::readDNAStringSet(fasta)
  gc_content <- compute_gc_content(sequence = sequence)
  
  gc_content <- reshape2::melt(gc_content)
  gc_content <- cbind.data.frame('sequence_id' = rownames(gc_content)
                                 ,  'gc_content'  = gc_content[,'value'])
  
  if( !is.null(outfile) ) {
    if( grepl('\\.rds', outfile) ) {
      saveRDS(gc_content, file = outfile)
    } else {
      save(gc_content, file = outfile)
    }
  }
  
  return(gc_content)
}

# CpG content ====
compute_expected_CpG <- function(sequence
                                 , method = 'G-G-F' # Gardiner-Frommer / (S-B-B) Saxonov-Berg-Brutlag
)
{
  n <- sapply(Biostrings::DNA_BASES, Biostrings::vcountPattern, subject=sequence)
  rownames(n) <- names(sequence)
  
  if( method == 'G-G-F' ) {
    expected_cpg <- ( n[,'C']*n[,'G'] )/rowSums(n)
    
  } else if( method == 'S-B-B') {
    expected_cpg <- ((( n[,'G'] + n[,'C'] )/2)**2)/rowSums(n)
    # expected_cpg <- (compute_gc_content(sequence = sequence )/2)**2
    
  } else {
    stop(message('[!] Please, provide correct method:'
                 ,'\n --- (G-G-F) Gardiner-Garden-Frommer'
                 ,'\n --- (S-B-B) Saxonov-Berg-Brutlag'))
  }
  
  return(expected_cpg)
}

compute_normalized_CpG <- function(sequence, ...)
{
  
  n <- sapply(Biostrings::DNA_BASES, Biostrings::vcountPattern, subject=sequence)
  
  observed_cpg <- Biostrings::vcountPattern(pattern = 'CG', subject = sequence)
  names(observed_cpg) <- names(sequence)
  expected_cpg <- compute_expected_CpG(sequence = sequence, ...)
  normalized_cpg <- observed_cpg/expected_cpg
  
  return(normalized_cpg)
  
}

get_normCpG_from_fasta <- function(fasta
                                   , outfile = NULL
                                   , ...)
{
  
  message("[+] Compute normalized CpG content from fasta")
  message("  -- file: "   , fasta)
  
  sequence   <- Biostrings::readDNAStringSet(fasta)
  normCpG_content <- compute_normalized_CpG(sequence = sequence, ...)
  
  normCpG_content <- reshape2::melt(normCpG_content)
  normCpG_content <- cbind.data.frame('sequence_id' = rownames(normCpG_content)
                                      ,  'CpG o/e ratio'  = normCpG_content[,'value'])
  
  if( !is.null(outfile) ) {
    if( grepl('\\.rds', outfile) ) {
      saveRDS(normCpG_content, file = outfile)
    } else {
      save(normCpG_content, file = outfile)
    }
  }
  
  return(normCpG_content)
}

plot_cpg <- function(counts
                     , ptitle = NULL)
{
  
  counts <- reshape2::melt(counts, id.var='sequence_id')
  n <- length(unique(counts$variable))
  
  bw <- function(x) (2 * IQR(x) / length(x)^(1/3)) # Freedman–Diaconis rule
  
  p <- ggplot(counts, aes(x=value, fill=variable)) + facet_wrap(~variable, ncol=unique(n), scales = 'free_x') +
    geom_histogram(col='black', size = 0.25, alpha = 0.8, binwidth = bw) +
    theme_bw() + my_theme + xlab(NULL)  + scale_fill_aaas()
  
  if(!is.null(ptitle)) {
    p <- p + ggtitle(ptitle) +
      theme(plot.title = element_text(size = 8, hjust = 0.5, face = 'bold'))
  }
  
  return(p)
}

