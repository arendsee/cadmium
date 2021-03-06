# ============================================================================
# Statistics
# ============================================================================

#' Gumbel density function
#'
#' @param x number
#' @param mu mean
#' @param s standard deviation
#' @export
dgumbel <- function(x, mu, s){
  z <- (mu - x) / s
  exp( z-exp(z) ) / s
}

#' Gumbel function
#'
#' @param q quantile 
#' @param mu mean
#' @param s standard deviation
#' @export
pgumbel <- function(q, mu, s){
  z <- (q - mu) / s
  exp(-exp(-z))
}

#' Gumbel function
#'
#' @param p quantile 
#' @param mu mean
#' @param s standard deviation
#' @export
qgumbel <- function(p, mu, s){
  mu - s*log(-log(p))
}


# Given simulated data, generate functions for assigning significance to
# scores.
#
# @param sam data.frame with columns 'query', 'score', and 'logmn'
fit.gumbel <- function(sam){

  stopifnot(c('query', 'score', 'logmn') %in% names(sam))
  sam <- sam                         %>%
    # Filter out maximum score for each query
    dplyr::group_by(query)           %>%
    dplyr::filter(score==max(score)) %>%
    dplyr::summarize(score=mean(score), logmn=max(logmn)) # to remove ties

  fit <- L1pack::l1fit(sam$logmn, sam$score)$coefficients
  b0 <- fit[[1]]
  b1 <- fit[[2]]

  # # adj.fit <- robustreg::robustRegBS(formula = score ~ logmn, data = sam)
  # adj.fit <- robustreg::robustRegH(formula = score ~ logmn, data = sam)
  #
  # b0 <- coef(adj.fit)[1]
  # b1 <- coef(adj.fit)[2]

  get.adj.from.score <- function(x, logmn){
    x - b0 - b1 * logmn
  }

  get.score.from.adj <- function(x, logmn){
    x + b0 + b1 * logmn
  }

  scores <- get.adj.from.score(sam$score, sam$logmn)

  gumbel.fit <- fitdistrplus::fitdist(
    data   = scores,
    distr  = "gumbel",
    start  = list(mu=mean(scores), s=sd(scores)),
    method = "mge",
    gof    = "CvM" # mge distance method 
  )

  mu <- gumbel.fit$estimate['mu']
  s  <- gumbel.fit$estimate['s']

  p <- function(q, logmn){
    q <- get.adj.from.score(q, logmn)
    # from the Gumbel CDF: e^{-e^{-(q-mu)/s)}}
    z <- (q - mu) / s
    exp(-exp(-z))
  }

  q <- function(p, logmn){
    adj.score <- mu - s*log(-log(p))
    get.score.from.adj(adj.score, logmn)
  }

  d <- function(x, logmn){
    x <- get.adj.from.score(x, logmn)
    z <- (mu - x) / s
    exp( z-exp(z) ) / s
  }

  list(fit=gumbel.fit, p=p, d=d, q=q)
}


# ============================================================================
# Search sequence - AA-AA - Query genes against target genes
# ============================================================================

nothing <- function(x) NULL

aln_xy <- function(x, y, substitution_matrix="BLOSUM80", simulation=FALSE){

  aln <- Biostrings::pairwiseAlignment(
    pattern=x,
    subject=y,
    type='local',
    substitutionMatrix=substitution_matrix
  )
  S4Vectors::metadata(aln)$query  <- names(x) 
  S4Vectors::metadata(aln)$target <- names(y) 

  a <- data.frame(
    query  = names(x),
    target = names(y),
    score  = BiocGenerics::score(aln),
    qwidth = Biostrings::width(x),
    twidth = Biostrings::width(y),
    stringsAsFactors = FALSE
  )

  map <- dplyr::group_by(a, query) %>%
    # Calculate adjusted score
    dplyr::summarize(logmn=log2(qwidth[1]) + log2(sum(twidth))) %>%
    base::merge(a) %>%
    dplyr::select(query, target, score, logmn)

  list(map=map, aln=aln)
}

AA_aln <- function(queseq, tarseq, nsims=10000, ...){

  # Store the original query to target mapping for later testing.
  # The input and output must have the same mapping.
  original.pairs = data.frame(
    query=names(queseq),
    target=names(tarseq),
    stringsAsFactors=FALSE
  )
  original.length = length(queseq)

  aln_result <- aln_xy(queseq, tarseq, ...)
  map <- aln_result$map

  # Simulate best hit for each query against randomized and reversed target sequences
  # 1. sample number of target sequences per query
  times <- names(queseq)  %>%
      factor              %>%
      summary(maxsum=Inf) %>%
      as.numeric          %>%
      sample(nsims, replace=TRUE) 
  # 2. randomly select query ids
  simids <- sample(1:length(queseq), nsims, replace=TRUE) %>%
    # 3. use each a number of times drawn from the empirical targets per query
    # distribution
    rep(times=times)
  # 4. give each simulated query a unique id
  simnames <- paste0('t', 1:nsims) %>% rep(times=times)
  # 5. align these against random targets
  simulation_result <- aln_xy(
    queseq[simids] %>% set_names(simnames),
    tarseq %>%
        sample(length(simnames), replace=TRUE) %>%
        IRanges::reverse(),
    simulation=TRUE,
    ...
  )

  sam <- simulation_result$map %>% dplyr::select(-target)

  gum <- fit.gumbel(sam)

  map$pval <- 1 - gum$p(map$score, map$logmn)
  sam$pval <- 1 - gum$p(sam$score, sam$logmn)

  # Adjust returned sample size to size of query
  if(nlevels(sam$query) > nlevels(map$query)){
    simnames <- levels(sam$query) %>% sample(nlevels(map$query))
    sam <- subset(sam, query %in% simnames) %>% droplevels
  }

  # Assert the query to target mapping has not changed
  stopifnot(
    map            %>% dplyr::arrange(query, target) %$% target ==
    original.pairs %>% dplyr::arrange(query, target) %$% target
  )
  # Assert table length has not changed
  stopifnot(nrow(map) == original.length)

  list(
    map     = map,
    dis     = gum,
    sam     = sam,
    nsims   = nsims,
    aln     = aln_result$aln
  )
}

align_by_map <- function(
  queseq,
  tarseq,
  map,
  queries = names(queseq),
  permute = FALSE,
  ...
){

  "Queries may be missing from the map if there are no target genes in any of
  their search intervals"

  map <- dplyr::select(map, .data$query, .data$target)

  if(permute){
    map$query <- map$query %>% { .[sample.int(length(.))] }
  }

  map <- map[map$query %in% queries, ]

  if(!all(map$target %in% names(tarseq))){
    dif <- setdiff(map$target, names(tarseq))
    msg <- "%s of %s of target genes missing in map; [%s, ...]"
    warning(sprintf(
      msg,
      length(dif),
      length(tarseq),
      paste(head(dif), collapse=', ')
    ))

    map <- map[map$target %in% names(tarseq), ]
    map <- map[map$query %in% names(queseq),  ]
  }

  # Pair the focal gene to each target sequence
  queseq <- queseq[ match(map$query,  names(queseq)) ]
  tarseq <- tarseq[ match(map$target, names(tarseq)) ]

  AA_aln(queseq, tarseq, ...)
}



# ============================================================================
# DNA alignment
# ============================================================================

get_logmn <- function(pattern, subject){
  data.frame(
    seqid = names(pattern),
    qlen = BiocGenerics::width(pattern),
    tlen = BiocGenerics::width(subject)
  ) %>%
    dplyr::group_by(seqid) %>%
    dplyr::mutate(logmn = log2(qlen[1]) + log2(sum(tlen))) %$% logmn 
}

alignToGenome <- function(
  queseq,
  tarseq,
  simulation = FALSE,
  offset     = 0,
  permute    = FALSE
){

    # Search + and - strands
    pattern <- c(queseq, Biostrings::reverseComplement(queseq))
    subject <- c(tarseq, tarseq)

    # If this is a control, permute the indices
    # This is used to determine hit significance
    if(permute){
      permid <- sample.int(length(subject))
      subject <- subject[permid]
      # Reverse the base order (but DO NOT complement)
      subject <- IRanges::reverse(subject)
      if(length(offset) == length(tarseq)){
        offset <- c(offset, offset)[permid]
      }
    }

    aln <- Biostrings::pairwiseAlignment(
      pattern=pattern,
      subject=subject,
      type='local'
    )
    S4Vectors::metadata(aln)$query  <- names(queseq) 
    S4Vectors::metadata(aln)$target <- names(tarseq) 

    CNEr::GRangePairs(
      first = GenomicRanges::GRanges(
        seqnames = names(pattern),
        ranges = IRanges::IRanges(
          start = Biostrings::pattern(aln) %>% BiocGenerics::start(),
          end = Biostrings::pattern(aln) %>% BiocGenerics::end()
        )
      ),
      second = GenomicRanges::GRanges(
        seqnames = names(subject),
        ranges = IRanges::IRanges(
          start = Biostrings::subject(aln) %>% BiocGenerics::start() %>% '+'(offset),
          end = Biostrings::subject(aln) %>% BiocGenerics::end() %>% '+'(offset)
        )
      ),
      strand = c(rep('+', length(aln)/2), rep('-', length(aln)/2)),
      score = BiocGenerics::score(aln),
      logmn = get_logmn(pattern, subject),
      query = names(pattern),
      qwidth = Biostrings::width(pattern),
      twidth = Biostrings::width(subject)
    )

}


get_dna2dna <- function(queseq, tarseq, queries, offset, maxspace=1e8){

  "Currently sequence searching is performed with a quadratic time alignment
  algorithm. If the search space is too big, this becomes prohibitively
  time-consuming. Ultimately I will need to use a heuristic program, such as
  BLAST. But this is not yet implemented. For now I just ignore spaces that
  are too large. Specifically, if `n*m > maxspace`, I skip the pair. The
  skipped pairs will be returned and ultimately be represented as a column in
  the feature table."

  # The query sequences are already paired against the search sequences
  stopifnot(length(queseq) == length(tarseq))
  # There exists an offset for each search interval
  stopifnot(length(tarseq) == length(offset))
  # An offset of 0 indicates the search interval starts at start of scaffold
  stopifnot(offset >= 0)
  # All query names are in the full list of focal genes
  stopifnot(queries %in% names(queseq))

  too_big <- log(GenomicRanges::width(queseq)) +
             log(GenomicRanges::width(tarseq)) > log(maxspace)

  skipped <- names(queseq)[too_big] %>% unique

  i <- setdiff(queries, skipped) %>% base::match(names(queseq))
  queseq <- queseq[i]
  tarseq <- tarseq[i]
  offset <- offset[i]

  ctrl <- alignToGenome(
            queseq     = queseq,
            tarseq     = tarseq,
            offset     = offset,
            permute    = TRUE,
            simulation = TRUE
          )

  hits <- alignToGenome(
            queseq     = queseq,
            tarseq     = tarseq,
            offset     = offset,
            permute    = FALSE,
            simulation = FALSE
          )

  gum <-
    S4Vectors::mcols(ctrl) %>%
    as.data.frame %>%
    dplyr::select(.data$query, .data$score, .data$logmn) %>%
    fit.gumbel

  S4Vectors::mcols(hits)$pval <- 1 - gum$p(S4Vectors::mcols(hits)$score,
                                           S4Vectors::mcols(hits)$logmn)
  S4Vectors::mcols(ctrl)$pval <- 1 - gum$p(S4Vectors::mcols(ctrl)$score,
                                           S4Vectors::mcols(ctrl)$logmn)

  list(
    map     = hits,
    sam     = ctrl,
    dis     = gum,
    skipped = skipped
  )
}
