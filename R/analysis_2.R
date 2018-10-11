# Collect secondary data for one species. The main output is a feature table.
# Internal data is summarized (for use in downstream diagnostics or
# side-analyses) and cached.
compare_target_to_focal <- function(m, con, species, group, gene_tag){

  #- Overlaps :: Table {
  #-    query  = Seqid
  #-  , target = Seqid
  #-  , type   = SOType
  #-  , qid    = ? -- I think this is an integer ID, but need to double check
  #-  , tid    = ?
  #-  }

  # TODO:
  # [ ] look into deeper synder analysis
  # [ ] formalize results in classes

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

  # SyntenyMap -> mRNA -> SearchIntervals
  .view_target(m, "synmap") %>%
    # fsi_
    rmonad::funnel(
      syn = .,
      gff = .view_focal(., "mRNA")
    ) %*>%
    synder::search(
      swap    = FALSE,
      trans   = con@synder@trans,
      k       = con@synder@k,
      r       = con@synder@r,
      offsets = con@synder@offsets
    ) %>% .tag("synder_out") %>>%
    # synder_summary_
    synder::flag_summary() %>% .tag("summary_synder_flags") %>%

  .view("synder_out") %>>%
    find_indels(indel.threshold=con@alignment@indel_threshold) %>%
    .tag("indels") %>%

  .view_target("nstring") %>>%
    synder::as_gff %>%
    rmonad::funnel(si = .view(., "synder_out"), gff=.) %*>% overlapMap %>%
    rmonad::funnel(
      x = .,
      nstrings = .view_target(., "nstring")
    ) %*>% {

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
  rmonad::funnel(
    si  = .view(., "synder_out"),
    gff = .view_target(., "mRNA")
  ) %*>%
    overlapMap %>% .tag("f_si_map") %>%

  # f_si_map_orf_
  .view_target("orfgff") %>>%
    # TODO: do this in analysis_1
    { synder::as_gff(., id=GenomicRanges::mcols(.)$seqid) } %>%
    rmonad::funnel(
      si  = .view(., "synder_out"),
      gff = .
    ) %*>% overlapMap %>% .tag("f_si_map_orf") %>%

  # - align query protein against target genes in the SI
  {
    rmonad::funnel(
      queseq  = .view_focal(., "faa"),
      tarseq  = .view_target(., "faa"),
      map     = .view(., "f_si_map"),
      queries = rmonad::view(., gene_tag)
    )
  } %*>% align_by_map %>% .tag("aa2aa") %>%

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
  } %*>% align_by_map %>% .tag("aa2orf") %>%
    # TODO: Resume from here:
    # [ ] see summary_orfaa, the minus-strand was mis-translated
    # [ ] fix names in tarseq, currently they are the scaffold names,
    #     they should be the ORF names (e.g. orf_1)

  # # rand_aa2orf_ <- rmonad::funnel(
  # #   queseq  = f_faa_,
  # #   tarseq  = t_orffaa_,
  # #   map     = f_si_map_orf_,
  # #   queries = seqids
  # # ) %*>% align_by_map(permute=TRUE)
  #
  # # transorf_map_
  # rmonad::funnel(
  #   transorfgff = .view_target(., "transorfgffDB"),
  #   f_si_map = .view(., "f_si_map")
  # ) %*>% {
  #   transorfgff %>%
  #   {
  #     data.frame(
  #       seqid = GenomicRanges::seqnames(.),
  #       orfid = GenomicRanges::mcols(.)$seqid,
  #       stringsAsFactors=FALSE
  #     )
  #   } %>%
  #   merge(f_si_map, by.x="seqid", by.y="target") %>%
  #   {
  #     data.frame(
  #       query = .$query,
  #       target = .$orfid,
  #       type = "transorf",
  #       stringsAsFactors=FALSE
  #     )
  #   }
  # } %>% .tag("f_si_map_transorf") %>%
  #
  # # aa2transorf - align query protein against ORFs in transcripts in the SI
  # rmonad::funnel(
  #   queseq  = .view_focal(., "faa"),
  #   tarseq  = .view_target(., "transorfaa"),
  #   map     = .view(., "f_si_map_transorf"),
  #   queries = rmonad::view(., gene_tag)
  # ) %*>% align_by_map %>% .tag('aa2transorf')
  #
  # # gene2genome - align query DNA sequence against the SI
  # rmonad::funnel(
  #   map     = .view(., "synder_out"),
  #   quedna  = .view_focal(., "transcript_seq"), # TODO: fix this in analysis_1
  #   tardna  = .view_target(., "genome_seq") ,
  #   queries = seqids
  # ) %*>% {
  #
  #   qids <- GenomicRanges::mcols(map)$attr
  #   if(! all(qids %in% names(quedna)) ){
  #
  #     ngood <- sum(unique(qids) %in% names(quedna))
  #     ntotal <- length(unique(qids))
  #
  #     msg <- "There is a mismatch between the names in the synteny map and
  #     those derived from the GFF file. %s of %s names from the synteny map are
  #     missing in the DNA file. Here are the first 5 names from the GFF-derived
  #     DNA file: [%s]. Here are the first 5 from the synteny map: [%s]."
  #
  #     stop(sprintf(
  #       msg,
  #       ngood,
  #       ntotal,
  #       paste0(head(names(quedna), 5), collapse=", "),
  #       paste0(head(qids, 5), collapse=", ")
  #     ))
  #
  #   }
  #
  #   tarseq <- Rsamtools::getSeq(x=tardna, CNEr::second(map))
  #   queseq <- quedna[ qids ]
  #   offset <- GenomicRanges::start(CNEr::second(map)) - 1
  #
  #   # An rmonad bound list with elements: map | sam | dis | skipped
  #   get_dna2dna(
  #     tarseq  = tarseq,
  #     queseq  = queseq,
  #     queries = queries,
  #     offset  = offset
  #   )
  #
  # } %>% .tag('gene2genome') %*>% rmonad::funnel(
  #   cds  = tcds_,
  #   exon = texons_,
  #   mrna = tgff_
  # ) %*>% {
  #
  #   "Find the queries that overlap a CDS, exon, or mRNA (technically pre-mRNA,
  #   but GFF uses mRNA to specify the entire transcript before splicing)"
  #
  #   met <- GenomicRanges::mcols(map)
  #   rng <- CNEr::second(map)
  #
  #   has_cds_match  <- GenomicRanges::findOverlaps(rng, cds)  %>% queryHits %>% unique
  #   has_exon_match <- GenomicRanges::findOverlaps(rng, exon) %>% queryHits %>% unique
  #   has_mrna_match <- GenomicRanges::findOverlaps(rng, mrna) %>% queryHits %>% unique
  #
  #   met$cds_match  <- seq_along(met$score) %in% has_cds_match
  #   met$exon_match <- seq_along(met$score) %in% has_exon_match
  #   met$mrna_match <- seq_along(met$score) %in% has_mrna_match
  #
  #   GenomicRanges::mcols(map) <- met
  #
  #   list(
  #     map=map,
  #     sam=sam,
  #     dis=dis,
  #     skipped=skipped
  #   )
  # }

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
      species  = target_species,
      group    = "control",
      gene_tag = "control_genes"
    )
  }

  m
}
