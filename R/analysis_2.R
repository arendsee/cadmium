get_gff_from_txdb <- function(x, ...){
  synder::as_gff(
    x,
    id=GenomicRanges::mcols(x)$tx_name,
    ...
  )
}

compare_target_to_focal <- function(
  species,
  fgff,
  f_primary,
  t_primary,
  synmap,
  seqids,
  con
){

  "
  Collect secondary data for one species. The main output is a feature table.
  Internal data is summarized (for use in downstream diagnostics or
  side-analyses) and cached.
  "

  ttxdb_ <- rmonad::as_monad(from_cache(t_primary@files@gff.file, type='sqlite'))

  tgff_ <- ttxdb_ %>>% {
    get_gff_from_txdb(
      x        = GenomicFeatures::transcripts(.),
      seqinfo_ = t_primary@seqinfo,
      type     = 'mRNA'
    )
  }

  tcds_ <- ttxdb_ %>>% {
    get_gff_from_txdb(
      x        = GenomicFeatures::cds(.),
      seqinfo_ = t_primary@seqinfo,
      type     = 'CDS'
    )
  }

  texons_ <- ttxdb_ %>>% {
    get_gff_from_txdb(
      x        = GenomicFeatures::exons(.),
      seqinfo_ = t_primary@seqinfo,
      type     = 'exon'
    )
  }

  synmap_ <- rmonad::as_monad( from_cache(synmap@synmap.file) )

  fsi_ <- rmonad::funnel(
    syn = synmap_,
    gff = fgff
  ) %*>%
  synder::search(
    swap    = FALSE,
    trans   = con@synder@trans,
    k       = con@synder@k,
    r       = con@synder@r,
    offsets = con@synder@offsets
  )

  ## TODO: resurrect?
  # esi_ <- rmonad::funnel(
  #   syn = synmap_,
  #   gff = fsi_ %>>% {
  #
  #     "Echo the search intervals back at the query."
  #
  #     synder::as_gff(
  #       CNEr::second(.),
  #       id=as.character(seq_along(.))
  #     )
  #   }
  # ) %*>%
  # synder::search(
  #   swap    = TRUE,
  #   trans   = con@synder@trans,
  #   k       = con@synder@k,
  #   r       = con@synder@r,
  #   offsets = con@synder@offsets
  # )

  ## TODO: resurrect?
  # rsi_ <- rmonad::funnel(
  #   syn = synmap_,
  #   gff = tgff_
  # ) %*>% {
  #
  # "
  # Trace target genes to focal genome using the reversed synteny map. This may
  # produce 'SKIPPING ENTRY' warnings. These warnings mean that the GFF
  # contains scaffolds that are missing in the synteny map. This is likely to
  # occur in genomes that are not fully assembled (thousands of scaffolds).
  # "
  #
  #   synder::search(
  #     syn,
  #     gff,
  #     swap    = TRUE,
  #     trans   = con@synder@trans,
  #     k       = con@synder@k,
  #     r       = con@synder@r,
  #     offsets = con@synder@offsets
  #   )
  # }

  synder_summary_ <- fsi_ %>>% synder::flag_summary

  indels_ <- rmonad::funnel(fsi_, con@alignment@indel_threshold) %*>% find_indels

  nstrings_ <- from_cache(t_primary@files@nstring.file) %>>% synder::as_gff

  gapped_ <- rmonad::funnel(
      x = rmonad::funnel(si = fsi_, gff = nstrings_) %*>% overlapMap
    , nstrings = nstrings_
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
  }

  # TODO: there is a lot more analysis that could be done here ...

  # f_count_qid_ <- f_si_map_ %>>% {
  #   "The number of target genes in each search interval"
  #   dplyr::group_by(., .data$qid) %>% dplyr::count()
  # }
  #
  # f_count_tid_ <- f_si_map_ %>>% {
  #   "The number of search intervals that overlap given target gene"
  #   dplyr::group_by(., .data$tid) %>% dplyr::count()
  # }

  f_si_map_ <- rmonad::funnel(si=fsi_, gff=tgff_) %*>% overlapMap

  # r_si_map_ <- rmonad::funnel(si=rsi_, gff=fgff) %*>% overlapMap

  orfgff_ <- t_primary@files@orfgff.file %>>%
    from_cache %>>%
    { synder::as_gff(., id=GenomicRanges::mcols(.)$seqid) }

  f_si_map_orf_ <- rmonad::funnel(
    si  = fsi_,
    gff = orfgff_
  ) %*>% overlapMap

  f_faa_ <- f_primary@files@aa.file %>>% from_cache
  t_faa_ <- t_primary@files@aa.file %>>% from_cache

  t_orffaa_ <- rmonad::as_monad(t_primary@files@orffaa.file %>% from_cache)

  t_transorffaa_ <- rmonad::as_monad(t_primary@files@transorffaa.file %>% from_cache)

  # - align query protein against target genes in the SI
  aa2aa_ <- rmonad::funnel(
    queseq  = f_faa_,
    tarseq  = t_faa_,
    map     = f_si_map_,
    queries = seqids
  ) %*>% align_by_map

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
  aa2orf_ <- rmonad::funnel(
    queseq  = f_faa_,
    tarseq  = t_orffaa_,
    map     = f_si_map_orf_,
    queries = seqids
  ) %*>% align_by_map

  # rand_aa2orf_ <- rmonad::funnel(
  #   queseq  = f_faa_,
  #   tarseq  = t_orffaa_,
  #   map     = f_si_map_orf_,
  #   queries = seqids
  # ) %*>% align_by_map(permute=TRUE)

  transorf_map_ <- rmonad::funnel(
    transorfgff = t_primary@files@transorfgff.file %>% from_cache,
    f_si_map = f_si_map_
  ) %*>% {
    transorfgff %>%
    {
      data.frame(
        seqid = GenomicRanges::seqnames(.),
        orfid = GenomicRanges::mcols(.)$seqid,
        stringsAsFactors=FALSE
      )
    } %>%
    merge(f_si_map, by.x='seqid', by.y='target') %>%
    {
      data.frame(
        query = .$query,
        target = .$orfid,
        type = 'transorf',
        stringsAsFactors=FALSE
      )
    }
  }

  # - align query protein against ORFs in transcripts in the SI
  aa2transorf_ <- rmonad::funnel(
    queseq  = f_faa_,
    tarseq  = t_transorffaa_,
    map     = transorf_map_,
    queries = seqids
  ) %*>% align_by_map


  # - align query DNA sequence against the SI
  gene2genome_ <- rmonad::funnel(
    map     = fsi_,
    quedna  = f_primary@files@trans.file %>>% from_cache %>>% Rsamtools::scanFa,
    tardna  = t_primary@files@dna.file %>>% from_cache,
    queries = seqids
  ) %*>% {

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

    # An rmonad bound list with elements: map | sam | dis | skipped
    get_dna2dna(tarseq=tarseq, queseq=queseq, queries=queries, offset=offset)

  } %*>% rmonad::funnel(
    cds  = tcds_,
    exon = texons_,
    mrna = tgff_
  ) %*>% {

    "Find the queries that overlap a CDS, exon, or mRNA (technically pre-mRNA,
    but GFF uses mRNA to specify the entire transcript before splicing)"

    met <- GenomicRanges::mcols(map)
    rng <- CNEr::second(map)

    has_cds_match  <- GenomicRanges::findOverlaps(rng, cds)  %>% queryHits %>% unique
    has_exon_match <- GenomicRanges::findOverlaps(rng, exon) %>% queryHits %>% unique
    has_mrna_match <- GenomicRanges::findOverlaps(rng, mrna) %>% queryHits %>% unique

    met$cds_match  <- seq_along(met$score) %in% has_cds_match
    met$exon_match <- seq_along(met$score) %in% has_exon_match
    met$mrna_match <- seq_along(met$score) %in% has_mrna_match

    GenomicRanges::mcols(map) <- met

    list(
      map=map,
      sam=sam,
      dis=dis,
      skipped=skipped
    )
  }

  rmonad::funnel(
    queries        = seqids,
    si             = fsi_,
    synder_summary = synder_summary_,
    indels         = indels_,
    # f_si_map       = f_si_map_,
    # r_si_map       = r_si_map_,
    aa2aa          = aa2aa_,
    aa2orf         = aa2orf_,
    aa2transorf    = aa2transorf_,
    gene2genome    = gene2genome_,
    gapped         = gapped_
  )

}

#' Load secondary data
#'
#' @export
secondary_data <- function(primary_input, con){

  "
  Loop over all target species, produce the secondary data needed to classify
  orphans. This step will take a long time (many CPU hours).
  "

  focal_species <- con@input@focal_species
  f_primary <- primary_input@species[[focal_species]]

  focal_gff_ <- from_cache(f_primary@files@gff.file, type='sqlite') %>>%
    GenomicFeatures::transcripts %>>%
    get_gff_from_txdb(seqinfo_=f_primary@seqinfo, type='mRNA') %>_%
    {

      "Assert query and control ids match the names in the GFF file"

      check_ids <- function(ids, label){
        txnames <- GenomicRanges::mcols(.)$attr
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

      check_ids(primary_input@queries, "query")
      check_ids(primary_input@control, "control")

    }

  rmonad::funnel(
    prim = primary_input,
    fgff = focal_gff_
  ) %*>% {

    focal_faa <- f_primary@files@aa.file %>>% from_cache %>% m_value

    die_if_genes_are_missing <- function(whole, part, label){
      missing <- !(part %in% whole)
      if(sum(missing) > 0){
        msg <- "%s of %s %s gene names are missing in the focal species (%s). These may be genes that are not translated (e.g. tRNAs or partial genes with no CDS). To diagnose, check the focal species GFF for the following missing ids: %s"
        stop(sprintf(
          msg,
          sum(missing),
          length(part),
          label,
          focal_species,
          paste(part[missing], collapse=", ")
        ))
      }
    }

    die_if_genes_are_missing(names(focal_faa), prim@queries, "query")
    die_if_genes_are_missing(names(focal_faa), prim@control, "control")

    species <- names(prim@species)
    target_species <- setdiff(species, focal_species)

    ss <- target_species %>%
      lapply(
        function(spec){
          rmonad::funnel(
            query_results = compare_target_to_focal(
              species   = spec,
              fgff      = fgff,
              f_primary = f_primary,
              t_primary = prim@species[[spec]],
              synmap    = prim@synmaps[[spec]],
              seqids    = prim@queries,
              con       = con
            ),
            control_results = compare_target_to_focal(
              species   = spec,
              fgff      = fgff,
              f_primary = f_primary,
              t_primary = prim@species[[spec]],
              synmap    = prim@synmaps[[spec]],
              seqids    = prim@control,
              con       = con
            )
          )
        }
      )

    names(ss) <- target_species
    rmonad::combine(ss)

  }

}
