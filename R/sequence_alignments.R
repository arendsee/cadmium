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
  adj.fit <- robustreg::robustRegBS(score ~ logmn, sam)
  b0 <- coef(adj.fit)[1]
  b1 <- coef(adj.fit)[2]

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
    method = "mge",gof="CvM" # mge distance method 
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

aln_xy <- function(x, y){
  a <- data.frame(
    query = names(x),
    target = names(y),
    score = Biostrings::pairwiseAlignment(
      pattern=x,
      subject=y,
      type='local',
      substitutionMatrix='BLOSUM80',
      scoreOnly=TRUE
    ),
    qwidth=width(x),
    twidth=width(y),
    stringsAsFactors=FALSE
  )
  dplyr::group_by(a, query) %>%
    # Calculate adjusted score
    dplyr::summarize(logmn=log2(qwidth[1]) + log2(sum(twidth))) %>%
    base::merge(a) %>%
    dplyr::select(query, target, score, logmn)
}

AA_aln <- function(queseq, tarseq, nsims=10000){

  # Store the original query to target mapping for later testing.
  # The input and output must have the same mapping.
  original.pairs = data.frame(
    query=names(queseq),
    target=names(tarseq),
    stringsAsFactors=FALSE
  )
  original.length = length(queseq)

  map <- aln_xy(queseq, tarseq)

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
  sam <- aln_xy(
    queseq[simids] %>% set_names(simnames),
    tarseq %>% base::sample(length(simnames), replace=TRUE) %>% reverse
  ) %>% 
  dplyr::select(-target)

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
    map   = map,
    dis   = gum,
    sam   = sam,
    nsims = nsims
  )
}

align_by_map <- function(
  queseq,
  tarseq,
  map,
  queries = names(queseq),
  permute = FALSE
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

  AA_aln(queseq, tarseq, nsims=1000)

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

get_dna2dna <- function(queseq, tarseq, maxspace=1e8){

  set.seed(42)

  too.big <- log(width(queseq)) + log(width(tarseq)) > log(maxspace)
  if(any(too.big)){
    warning(sprintf('%d(%.1f%%) query / SI pairs are very large, N*M>%d. These
    pairs are ignored. Dealing with them will require a heuristic aligment
    program, such as BLAST, which is not currently implemented.',
    sum(too.big), signif(100*sum(too.big)/length(too.big), 2), maxspace))
    orfgff <- orfgff[!too.big]
    skipped <- ogen[too.big] %>% unique %>% names
    ogen <- ogen[!too.big]
  } else {
    skipped = c()
  }

  message(sprintf('This may require on the order of %.1f minutes',
    (wdith(queseq) * width(tarseq) * (3/9e8)) %>% sum %>% signif(1)))

  # Search + and - strands
  nuc.scores <- pairwiseAlignment(
    pattern=c(queseq, reverseComplement(queseq)),
    subject=c(tarseq, tarseq),
    type='local',
    scoreOnly=TRUE
  )

  hits <- data.frame(
    query            = queseq %>% names %>% rep(2),
    qwidth           = queseq %>% width %>% rep(2),
    twidth           = tarseq %>% width %>% rep(2),
    score            = nuc.scores,
    strand           = rep(c('+', '-'), each=length(query.seqs)),
    stringsAsFactors = FALSE
  ) %>%
  add_logmn 

  # Align queries against random search intervals

  ctrl <- pairwiseAlignment(
    pattern=c(queseq, reverseComplement(queseq)),
    subject=c(tarseq, tarseq),
    type='local',
    scoreOnly=TRUE
  )                 %>%
    add_logmn       %>%
    group_by(query) %>%
    filter(score == max(score))

  gum <- fit.gumbel(ctrl)

  hits$pval <- 1 - gum$p(hits$score, hits$logmn)
  ctrl$pval <- 1 - gum$p(ctrl$score, ctrl$logmn)

  list(
    map      = hits,
    dis      = gum,
    sam      = ctrl,
    skipped  = skipped,
    maxspace = maxspace
  )
}
