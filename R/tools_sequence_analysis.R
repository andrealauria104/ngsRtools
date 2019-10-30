# Tools for sequence analysis
require(Biostrings)
require(GenomicRanges)
# Utils ====
sequence_to_list <- function(fasta)
{
  fa <- readDNAStringSet(fasta)

  sequence_list <- vector(mode = 'list', length = length(fa))
  names(sequence_list) <- names(fa)

  for( i in names(fa)) {
    sequence_list[[i]] <- toString(fa[[i]])
  }

  rm(fa)
  gc(verbose = F)

  return(sequence_list)
}

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

# GC content ====
compute_gc_content <- function(sequence)
{
  require(Biostrings)

  n <- sapply(DNA_BASES, vcountPattern, subject=sequence)
  rownames(n) <- names(sequence)
  gc_content  <- round((n[,'G']+n[,'C'])/rowSums(n),4)

  return(gc_content)
}

get_gc_from_fasta <- function(fasta
                              , outfile = NULL
                              # , ...
                              )
{

  require(parallel)

  message("[+] Compute GC content from fasta")
  message("  -- file: "   , fasta)

  sequence   <- readDNAStringSet(fasta)
  gc_content <- compute_gc_content(sequence = sequence)
  # sequence_list  <- sequence_to_list(fasta = fasta)
  #
  # gc_content <- mclapply(sequence_list
  #                        , compute_gc_content
  #                        , mc.silent = T
  #                        , ... )
  #
  # gc_content <- reshape2::melt(gc_content)
  # colnames(gc_content) <- c('gc_content', 'sequence_id')

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
  require(Biostrings)

  n <- sapply(DNA_BASES, vcountPattern, subject=sequence)
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
  require(Biostrings)

  n <- sapply(DNA_BASES, vcountPattern, subject=sequence)

  observed_cpg <- vcountPattern(pattern = 'CG', subject = sequence)
  names(observed_cpg) <- names(sequence)
  expected_cpg <- compute_expected_CpG(sequence = sequence, ...)
  normalized_cpg <- observed_cpg/expected_cpg

  return(normalized_cpg)

}

get_normCpG_from_fasta <- function(fasta
                              , outfile = NULL
                              , ...)
{

  require(parallel)

  message("[+] Compute normalized CpG content from fasta")
  message("  -- file: "   , fasta)

  # sequence_list  <- sequence_to_list(fasta = fasta)
  #
  # normCpG_content <- mclapply(sequence_list
  #                        , compute_normalized_CpG
  #                        , mc.silent = T
  #                        , ... )
  sequence   <- readDNAStringSet(fasta)
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

  #source("theme_setting.R")
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


# GenomicRanges Tools ====
extend <- function(x, upstream=0, downstream=0)
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}
