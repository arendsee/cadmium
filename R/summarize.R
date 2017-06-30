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
    shapiro_norm     = stats::shapiro.test(sample(x,      size=1000, replace = TRUE)),
    log_shapiro_norm = stats::shapiro.test(sample(log(x), size=1000, replace = TRUE)),
    density          = stats::density(x, kernel="gaussian")
  )
}

#' @rdname fagin_summary
#' @export
summarize_faa <- function(x){
# TODO: fill in the stub values
  new(
    "faa_summary",
    class_composition = NA_real_,
    # inherited from seq_summary
    seqids          = names(x),
    sizes           = nchar(x),
    nseq            = length(x),
    comp            = NA_real_,
    comp_dist       = list(),
    n_non_canonical = NA_integer_,
    n_initial_start = NA_integer_,
    n_terminal_stop = NA_integer_,
    length_sum      = summarize_numeric(nchar(x))
  )
}

#' @rdname fagin_summary
#' @export
summarize_fna <- function(x){
# TODO: fill in the stub values
  new(
    "fna_summary",
    nd50             = NA_integer_,
    n_triple         = NA_integer_,
    codon_preference = list(),
    # inherited from seq_summary
    seqids          = names(x),
    sizes           = nchar(x),
    nseq            = length(x),
    comp            = NA_real_,
    comp_dist       = list(),
    n_non_canonical = NA_integer_,
    n_initial_start = NA_integer_,
    n_terminal_stop = NA_integer_,
    length_sum      = summarize_numeric(nchar(x))
  )
}

#' @rdname fagin_summary
#' @export
summarize_gff <- function(x){

  seqstats <- dplyr::summarize(
    dplyr::group_by(x, .data$seqid),
    min   = min(.data$start),
    max   = max(.data$stop),
    mRNAs = sum(.data$type == "mRNA")
  )

  lengths <- lapply(c("mRNA", "CDS", "exon"),
    function(s) {
      d <- dplyr::filter(x, .data$type == s) %>%
           dplyr::transmute(length = .data$stop - .data$start + 1)
      summarize_numeric(d[[1]])
    }
  )

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
  width <- x$stop - x$start + 1
  summarize_numeric(width)
}

#' @rdname fagin_summary
#' @export
summarize_syn <- function(x){
  qwidth <- x$qstop - x$qstart + 1
  twidth <- x$tstop - x$tstart + 1
  new(
    "synmap_summary",
    nrow = nrow(x),
    query_target_log2_ratio = summarize_numeric(log2 (qwidth / twidth) ),
    width_summary           = summarize_numeric(qwidth),
    score_summary           = summarize_numeric(x$score)
  )
}
