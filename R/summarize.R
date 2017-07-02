#' Summarize functions
#'
#' @param x An object to summarized
#' @name fagin_summary
NULL


#' @rdname fagin_summary
#' @export
summarize_numeric <- function(x){
  new(
    "numeric_summary",
    min              = min(x),
    q25              = stats::quantile(x, probs=0.25),
    median           = stats::median(x),
    q75              = stats::quantile(x, probs=0.75),
    max              = max(x),
    mean             = mean(x),
    sd               = stats::sd(x),
    n                = length(x),
    density          = stats::density(x, kernel="gaussian")
  )
}

#' @rdname fagin_summary
#' @export
summarize_faa <- function(x){
  new(
    "faa_summary",
    class_composition = NA_real_,
    # inherited from seq_summary
    seqids          = names(x) %>% unique,
    sizes           = Biostrings::width(x),
    nseq            = length(x),
    comp            = Biostrings::alphabetFrequency(x),
    n_initial_start = sapply(x, grepl, pattern="^M") %>% sum,
    n_terminal_stop = sapply(x, grepl, pattern="*$") %>% sum
  )
}

#' @rdname fagin_summary
#' @export
summarize_dna <- function(x){
  new(
    "dna_summary",
    n_triple         = (Biostrings::width(x) %% 3 == 0) %>% sum,
    # inherited from seq_summary
    seqids          = names(x) %>% unique,
    sizes           = Biostrings::width(x),
    nseq            = length(x),
    comp            = Biostrings::alphabetFrequency(x),
    n_initial_start = sapply(x, grepl, pattern="^ATG") %>% sum,
    n_terminal_stop = sapply(x, grepl, pattern="(TAA|TGA|TAG)$", perl=TRUE) %>% sum
  )
}

#' @rdname fagin_summary
#' @export
summarize_gff <- function(x){

  xdf <- data.frame(
    seqid = x@seqnames,
    stop  = GenomicRanges::end(x),
    start = GenomicRanges::start(x),
    type  = mcols(x)$type 
  )

  seqstats <- dplyr::summarize(
    dplyr::group_by(xdf, .data$seqid),
    min   = min(.data$start),
    max   = max(.data$stop),
    mRNAs = sum(.data$type == "mRNA")
  )

  features <- c("mRNA", "CDS", "exon")

  lengths <- lapply(features,
    function(s) {
      d <- dplyr::filter(xdf, .data$type == s) %>%
           dplyr::transmute(length = .data$stop - .data$start + 1)
      summarize_numeric(d[[1]])
    }
  )
  names(lengths) <- features

  new(
    "gff_summary",
    seqstats    = seqstats,
    mRNA_length = lengths$mRNA,
    CDS_length  = lengths$CDS,
    exon_length = lengths$exon
  )
}

#' @rdname fagin_summary
#' @export
summarize_nstring <- function(x){
  if(length(x) != 0 && nrow(x) > 3) {
    width <- x$stop - x$start + 1
    summarize_numeric(width)
  } else {
    # TODO: need a better means of dealing with N-less genomes
    summary(x)
  }
}

#' @rdname fagin_summary
#' @export
summarize_syn <- function(x){
  qwidth <- x$qstop - x$qstart + 1
  twidth <- x$tstop - x$tstart + 1
  new(
    "synmap_summary",
    nrow = nrow(x),
    width                   = summarize_numeric(qwidth),
    score                   = summarize_numeric(x$score),
    query_target_log2_ratio = summarize_numeric(log2 (qwidth / twidth) )
  )
}
