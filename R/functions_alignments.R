# ============================================================================
# Statistics
# ============================================================================

dgumbel <- function(x, mu, s){
  z <- (mu - x) / s
  exp( z-exp(z) ) / s
}

pgumbel <- function(q, mu, s){
  z <- (q - mu) / s
  exp(-exp(-z))
}

qgumbel <- function(p, mu, s){
  mu - s*log(-log(p))
}

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

aln_xy <- function(x, y, group, label, simulation=FALSE){

  aln <- Biostrings::pairwiseAlignment(
    pattern=x,
    subject=y,
    type='local',
    substitutionMatrix='BLOSUM80'
  )
  metadata(aln)$query  <- names(x) 
  metadata(aln)$target <- names(y) 

  if(simulation){
    group <- paste0(group, "-sim")
  }

  alnfile <- to_cache(aln, group=group, label=label)

  a <- data.frame(
    query  = names(x),
    target = names(y),
    score  = score(aln),
    qwidth = width(x),
    twidth = width(y),
    stringsAsFactors = FALSE
  )

  list(
      map = dplyr::group_by(a, query) %>%
        # Calculate adjusted score
        dplyr::summarize(logmn=log2(qwidth[1]) + log2(sum(twidth))) %>%
        base::merge(a) %>%
        dplyr::select(query, target, score, logmn)
    , alnfile = alnfile
  )

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
    tarseq %>% base::sample(length(simnames), replace=TRUE) %>% reverse,
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
    alnfile = aln_result$alnfile
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

  queseq <- queseq[ match(map$query,  names(queseq)) ]
  tarseq <- tarseq[ match(map$target, names(tarseq)) ]

  AA_aln(queseq, tarseq, nsims=1000, ...)

}



# ============================================================================
# DNA alignment
# ============================================================================

add_logmn <- function(d){
  dplyr::group_by(d, query)                                     %>%
    # Calculate adjusted score
    dplyr::filter(twidth > 1 & qwidth > 1)                      %>%
    dplyr::summarize(logmn=log2(qwidth[1]) + log2(sum(twidth))) %>%
    base::merge(d)
}



alignToGenome <- function(
  queseq,
  tarseq,
  group,
  label,
  simulation = FALSE,
  offset     = 0,
  permute    = FALSE
){
  # Search + and - strands
  pattern <- c(queseq, Biostrings::reverseComplement(queseq))
  subject <- c(tarseq, tarseq)

  if(permute){
    permid <- sample.int(length(subject))
    subject <- subject[permid]
    if(length(offset) == length(tarseq)){
      offset <- c(offset, offset)[permid]
    }
  }

  aln <- Biostrings::pairwiseAlignment(
    pattern=pattern,
    subject=subject,
    type='local'
  )
  metadata(aln)$query  <- names(queseq) 
  metadata(aln)$target <- names(tarseq) 

  if(simulation){
    group = paste0(group, "-sim")
  }
  alnfile <- to_cache(aln, group=group, label=label)

  map_ <- aln %>>% {
    CNEr::GRangePairs(
      first = GenomicRanges::GRanges(
        seqnames = names(pattern),
        ranges = IRanges::IRanges(
          start = Biostrings::pattern(.) %>% BiocGenerics::start(),
          end = Biostrings::pattern(.) %>% BiocGenerics::end()
        )
      ),
      second = GenomicRanges::GRanges(
        seqnames = names(subject),
        ranges = IRanges::IRanges(
          start = Biostrings::subject(.) %>% BiocGenerics::start() %>% '+'(offset),
          end = Biostrings::subject(.) %>% BiocGenerics::end() %>% '+'(offset)
        )
      ),
      strand = c(rep('+', length(.)/2), rep('-', length(.)/2)),
      score = BiocGenerics::score(.),
      query = names(pattern),
      qwidth = Biostrings::width(pattern),
      twidth = Biostrings::width(subject)
    )
  }

  rmonad::funnel(
      map = map_
    , alnfile = alnfile
  )

}


add_logmn <- function(d){
  GenomicRanges::mcols(d) <-
    as.data.frame(GenomicRanges::mcols(d)) %>%
    dplyr::group_by(.data$query) %>%
    dplyr::mutate(logmn=log2(.data$qwidth) + log2(sum(.data$twidth)))
  d
}


get_dna2dna <- function(queseq, tarseq, queries, offset, maxspace=1e8, ...){

  skipped_ <- rmonad::as_monad({

    "Currently sequence searching is performed with a quadratic time alignment
    algorithm. If the search space is too big, this becomes prohibitively
    time-consuming. Ultimately I will need to use a heuristic program, such as
    BLAST. But this is not yet implemented. For now I just ignore spaces that
    are too large. Specifically, if `n*m > maxspace`, I skip the pair. The
    skipped pairs will be returned and ultimately be represented as a column in
    the feature table."

    too_big <- log(GenomicRanges::width(queseq)) +
               log(GenomicRanges::width(tarseq)) > log(maxspace)
    names(queseq)[too_big] %>% unique

  })

  truncated_seqs_ <-
    rmonad::funnel(
      skipped = skipped_,
      queries = queries
    ) %*>% {

      "Select which gene/search_interval pairs to align"

      i <- setdiff(queries, skipped) %>% match(names(queseq))

      list(
        queseq=queseq[i],
        tarseq=tarseq[i],
        offset=offset[i]
      )
    }

  # TODO: continue from here: need to extract the cached file
  ctrl_ <- truncated_seqs_ %*>%
    alignToGenome(
      permute    = TRUE,
      simulation = TRUE,
      ...
    ) %>>% { .$map <- add_logmn(.$map); . }
    ## TODO: And what the flip is the following commented code? Why is it still
    ## here? I kept it for some reason ...
    # dplyr::group_by(query) %>%
    # dplyr::filter(.data$score == max(.data$score))

  hits_ <- truncated_seqs_ %*>%
    alignToGenome(
      permute = FALSE,
      ...
    ) %>>% { .$map <- add_logmn(.$map); . }

  gum_ <- ctrl_ %>>%
    {
      GenomicRanges::mcols(.$map) %>%
      as.data.frame %>%
      dplyr::select(.data$query, .data$score, .data$logmn)
    } %>>%
    fit.gumbel

  rmonad::funnel(
    ctrl = ctrl_,
    hits = hits_,
    gum  = gum_
  ) %*>% {

    GenomicRanges::mcols(hits$map)$pval <- 1 - gum$p(GenomicRanges::mcols(hits$map)$score,
                                                     GenomicRanges::mcols(hits$map)$logmn)
    GenomicRanges::mcols(ctrl$map)$pval <- 1 - gum$p(GenomicRanges::mcols(ctrl$map)$score,
                                                     GenomicRanges::mcols(ctrl$map)$logmn)

    list(
      map=hits$map,
      sam=ctrl$map,
      dis=gum,
      alnfile=hits$alnfile
    )
  } %*>%
  rmonad::funnel(
    skipped=skipped_
  )

}
