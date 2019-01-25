# Collect secondary data for one species. The main output is a feature table.
# Internal data is summarized (for use in downstream diagnostics or
# side-analyses) and cached.
compare_target_to_focal <- function(m, con, species, group, gene_tag){

  message("Processing ", species, " ", group)

  .view_target <- function(m, tag){
    rmonad::view(m, tag, species)
  }

  .view_focal <- function(m, tag){
    rmonad::view(m, tag, con@input@focal_species)
  }

  .view <- function(m, tag){
    rmonad::view(m, tag, species, group)
  }

  .tag <- function(m, tag){
    rmonad::tag(m, tag, species, group)
  }

  # ---------------------------
  # dirty hack to get around a strange bug where some species loose their
  # Seqinfo objects in the mRNA field. 
  qseqinfo <- rmonad::get_value(m, tag=c('genomeSeq', con@input@focal_species))[[1]] %>%
    make_seqinfo(con@input@focal_species)
  tseqinfo <- rmonad::get_value(m, tag=c('genomeSeq', species))[[1]] %>%
    make_seqinfo(species)
  # ---------------------------

  nsims_prot <- con@alignment@simulation@prot2prot
  nsims_allorf <- con@alignment@simulation@prot2allorf
  nsims_transorf <- con@alignment@simulation@prot2transorf
  subMat <- con@alignment@substitutionMatrix

  # SyntenyMap -> mRNA -> SearchIntervals
  .view_target(m, "synmap") %>%
    # fsi_
    {
      rmonad::funnel(
        syn = .,
        gff = .view_focal(., "mRNA")
      )
    } %*>%
    synder::search(
      swap    = FALSE,
      trans   = con@synder@trans,
      k       = con@synder@k,
      r       = con@synder@r,
      tcl     = tseqinfo,
      qcl     = qseqinfo,
      offsets = con@synder@offsets
    ) %>% .tag("synder_out") %>>%
    # synder_summary_
    synder::flag_summary() %>% .tag("summary_synder_flags") %>%

  .view("synder_out") %>>%
    find_indels(indel.threshold=con@alignment@indel_threshold) %>%
    .tag("indels") %>%

  .view_target("nstring") %>%
    {
      rmonad::funnel(
        gff = .,
        si  = .view(., "synder_out")
      )
    } %*>% overlapMap %>%
    {
      rmonad::funnel(
        x = .,
        nstrings = .view_target(., "nstring")
      )
    } %*>% {

      "
      Get a N-string table with the columns

      query  - query id
      siid   - search interval id
      length - length of N-string

      Filter out any queries that overlap no N-strings.

      If there are no N-strings, return an empty dataframe.
      "

      if(is.null(nstrings)){
        data.frame(
          query  = character(0),
          siid   = integer(0),
          length = integer(0)
        )
      } else {
        data.frame(
          query  = x$query,
          siid   = x$qid,
          length = GenomicRanges::width(nstrings)[x$qid]
        ) %>% { .[!is.na(.$length), ] }
      }
    } %>% .tag("gaps") %>%

  # SeachIntervals -> GeneModels -> Overlaps
  {
    rmonad::funnel(
      gff = .view_target(., "mRNA"),
      si  = .view(., "synder_out")
    )
  } %*>%
    overlapMap %>% .tag("f_si_map") %>%

  # f_si_map_orf_
  {
    rmonad::funnel(
      gff = .view_target(., "orfgff"),
      si  = .view(., "synder_out")
    )
  } %*>% overlapMap %>% .tag("f_si_map_orf") %>%

  # - align query protein against target genes in the SI
  {
    rmonad::funnel(
      queseq  = .view_focal(., "faa"),
      tarseq  = .view_target(., "faa"),
      map     = .view(., "f_si_map"),
      queries = rmonad::view(., gene_tag)
    )
  } %*>% align_by_map(nsims=nsims_prot, substitutionMatrix=subMat) %>%
         .tag("aa2aa") %>%

  # # Run the exact test above, but with the query indices scrambled. For a
  # # reasonably sized genome, there should be very few hits (I should do the
  # # math and make a test)
  # rand_aa2aa_ <- rmonad::funnel(
  #   queseq  = f_faa_,
  #   tarseq  = t_faa_,
  #   map     = f_si_map_,
  #   queries = seqids
  # ) %*>% align_by_map(permute=TRUE)

  # - align query protein against ORFs in genomic intervals in SI
  {
    rmonad::funnel(
      queseq  = .view_focal(., "faa"),
      tarseq  = .view_target(., "orffaa"),
      map     = .view(., "f_si_map_orf"),
      queries = rmonad::view(., gene_tag)
    )
  } %*>% align_by_map(nsims=nsims_allorf, substitutionMatrix=subMat) %>%
         .tag("aa2orf") %>%

  # rand_aa2orf_ <- rmonad::funnel(
  #   queseq  = f_faa_,
  #   tarseq  = t_orffaa_,
  #   map     = f_si_map_orf_,
  #   queries = seqids
  # ) %*>% align_by_map(permute=TRUE)

  # transorf_map_
  {
    rmonad::funnel(
      transorfgff = .view_target(., "transorfgff"),
      f_si_map = .view(., "f_si_map")
    )
  } %*>% {
    transorfgff %>%
    {
      data.frame(
        seqid = GenomicRanges::seqnames(.),
        orfid = GenomicRanges::mcols(.)$attr,
        stringsAsFactors=FALSE
      )
    } %>%
    merge(f_si_map, by.x="seqid", by.y="target") %>%
    {
      data.frame(
        query = .$query,
        target = .$orfid,
        type = "transorf",
        stringsAsFactors=FALSE
      )
    }
  } %>% .tag("f_si_map_transorf") %>%

  # aa2transorf - align query protein against ORFs in transcripts in the SI
  {
    rmonad::funnel(
      queseq  = .view_focal(., "faa"),
      tarseq  = .view_target(., "transorfaa"),
      map     = .view(., "f_si_map_transorf"),
      queries = rmonad::view(., gene_tag)
    )
  } %*>% align_by_map(nsims=nsims_transorf, substitutionMatrix=subMat) %>%
         .tag('aa2transorf') %>%

  # gene2genome - align query DNA sequence against the SI
  {
    rmonad::funnel(
      map     = .view(., "synder_out"),
      quedna  = .view_focal(., "transcriptomeSeq"),
      tardna  = .view_target(., "genomeDB"),
      queries = rmonad::view(., gene_tag)
    )
  } %*>% {

    qids <- GenomicRanges::mcols(map)$attr
    if(! all(qids %in% names(quedna)) ){

      ngood <- sum(unique(qids) %in% names(quedna))
      ntotal <- length(unique(qids))

      msg <- "There is a mismatch between the names in the synteny map and
      those derived from the GFF file. %s of %s names from the synteny map are
      missing in the DNA file. Here are the first 5 names from the GFF-derived
      DNA file: [%s]. Here are the first 5 from the synteny map: [%s]."

      stop(sprintf(
        msg,
        ngood,
        ntotal,
        paste0(head(names(quedna), 5), collapse=", "),
        paste0(head(qids, 5), collapse=", ")
      ))

    }

    tarseq <- Rsamtools::getSeq(x=tardna, CNEr::second(map))
    queseq <- quedna[ qids ]
    offset <- GenomicRanges::start(CNEr::second(map)) - 1

    get_dna2dna(
      queseq   = queseq,
      tarseq   = tarseq,
      queries  = queries,
      offset   = offset,
      maxspace = con@alignment@dna2dna_maxspace
    )
  } %>% .tag('gene2genome_aln') %>% {
    rmonad::funnel(
      gene2genome = .,
      cds  = .view_target(., "mRNA"),
      exon = .view_target(., "exon"),
      mrna = .view_target(., "CDS")
    )
  } %*>% {

    "Find the queries that overlap a CDS, exon, or mRNA (technically pre-mRNA,
    but GFF uses mRNA to specify the entire transcript before splicing)"

    map <- gene2genome$map
    met <- S4Vectors::mcols(map)
    rng <- CNEr::second(map)

    has_cds_match  <- GenomicRanges::findOverlaps(rng, cds)  %>% S4Vectors::queryHits() %>% unique
    has_exon_match <- GenomicRanges::findOverlaps(rng, exon) %>% S4Vectors::queryHits() %>% unique
    has_mrna_match <- GenomicRanges::findOverlaps(rng, mrna) %>% S4Vectors::queryHits() %>% unique

    met$cds_match  <- seq_along(met$score) %in% has_cds_match
    met$exon_match <- seq_along(met$score) %in% has_exon_match
    met$mrna_match <- seq_along(met$score) %in% has_mrna_match

    S4Vectors::mcols(map) <- met

    map
  } %>% .tag('gene2genome')
}

query_control_gene_check <- function(gff, qgenes, cgenes) {

  "Assert query and control ids match the names in the GFF file. Also confirm
  that query and control ids do not overlap"

  check_ids <- function(ids, label){
    txnames <- GenomicRanges::mcols(gff)$attr
    matches <- ids %in% txnames
    if(!all(matches)){
      msg <- "%s of %s %s ids are missing in the focal GFF (%s). Here is the head: [%s]"
      stop(sprintf(
        msg,
        sum(!matches),
        length(ids),
        label,
        focal_species,
        paste(head(ids[!matches]), collapse=", ")
      ))
    }
  }

  shared <- intersect(qgenes, cgenes)
  if(length(shared) > 0){
    warning(sprintf("%s genes are shared between query and control"))
  }

  check_ids(qgenes, "query")
  check_ids(cgenes, "control")
}

#' Load secondary data
#'
#' Loop over all target species, produce the secondary data needed to classify
#' orphans. This step will take a long time (many CPU hours).
#'
#' @param m Rmonad pipeline from stage 2
#' @export
secondary_data <- function(m, con){

  .view_focal <- function(m, tag){
    rmonad::view(m, tag, con@input@focal_species)
  }

  m <- .view_focal(m, "mRNA") %>%
    rmonad::funnel(
      gff = .,
      qgenes = rmonad::view(., "query_genes"),
      cgenes = rmonad::view(., "control_genes")
    ) %*>%
    query_control_gene_check

  for(target_species in get_targets(con)){
    m <- compare_target_to_focal(
      m,
      con,
      species  = target_species,
      group    = "query",
      gene_tag = "query_genes"
    )
    m <- compare_target_to_focal(
      m,
      con,
      species  = target_species,
      group    = "control",
      gene_tag = "control_genes"
    )
  }

  m
}
