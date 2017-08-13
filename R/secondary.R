get_transcripts_from_txdb <- function(spec){

  from_cache(spec@files@gff.file, type='sqlite') %>>%
  GenomicFeatures::transcripts %>>%
  {
    synder::as_gff(
      .,
      seqinfo=spec@seqinfo,
      id=GenomicRanges::mcols(.)$tx_name,
      type="mRNA"
    )
  }
}



compare_target_to_focal <- function(
  species,
  fgff,
  f_primary,
  t_primary,
  synmap,
  queries,
  con
){

  "
  Collect secondary data for one species. The main output is a feature table.
  Internal data is summarized (for use in downstream diagnostics or
  side-analyses) and cached.
  "

  tgff_ <- get_transcripts_from_txdb(t_primary)

  synmap_ <- as_monad( from_cache(synmap@synmap.file) )

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


  esi_ <- rmonad::funnel(
    syn = synmap_,
    gff = fsi_ %>>% {

      "Echo the search intervals back at the query."

      synder::as_gff(
        CNEr::second(.),
        id=as.character(seq_along(.))
      )
    }
  ) %*>%
  synder::search(
    swap    = TRUE,
    trans   = con@synder@trans,
    k       = con@synder@k,
    r       = con@synder@r,
    offsets = con@synder@offsets
  )

  rsi_ <- rmonad::funnel(
    syn = synmap_,
    gff = tgff_
  ) %*>% {

  "
  Trace target genes to focal genome using the reversed synteny map. This may
  produce 'SKIPPING ENTRY' warnings. These warnings mean that the GFF
  contains scaffolds that are missing in the synteny map. This is likely to
  occur in genomes that are not fully assembled (thousands of scaffolds).
  "

    synder::search(
      syn,
      gff,
      swap    = TRUE,
      trans   = con@synder@trans,
      k       = con@synder@k,
      r       = con@synder@r,
      offsets = con@synder@offsets
    )
  }

  synder_flags_summary_ <-
    rmonad::funnel(
      si      = fsi_,
      queries = queries
    ) %*>% summarize_syntenic_flags

  unassembled_ <- si_ %>>% find_unassembled

  scrambled_ <- si_ %>>% find_scrambled

  indels_ <- rmonad::funnel(si_, con@alignment@indel_threshold) %*>% find_indels

  nstrings_ <- from_cache(t_primary@files@nstring.file) %>>% synder::as_gff

  gapped_ <- rmonad::funnel(
    si = fsi_,
    gff = nstrings_
  ) %*>% overlapMap %>%
  { rmonad::funnel(x=., nstring=nstrings_) } %*>%
  {

    "
    Get a N-string table with the columns

    query  - query id
    siid   - search interval id
    length - length of N-string

    Filter out any queries that overlap no N-strings.

    If there are no N-strings, return an empty dataframe.
    "

    if(is.null(nstring)){
      data.frame(
        query  = character(0),
        siid   = integer(0),
        length = integer(0)
      )
    } else {
      data.frame(
        query  = x$query,
        siid   = x$qid,
        length = GenomicRanges::width(nstring)[x$qid]
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

  # rmonad::funnel(fmap=f_si_map_, rmap=r_si_map_)

  f_si_map_ <- rmonad::funnel(si=fsi_, gff=tgff_) %*>% overlapMap

  r_si_map_ <- rmonad::funnel(si=rsi_, gff=fgff) %*>% overlapMap

  orfgff_ <- t_primary@files@orfgff.file %>>%
    from_cache %>>%
    { synder::as_gff(., id=GenomicRanges::mcols(.)$seqid) }

  f_si_map_orf_ <- rmonad::funnel(
    si  = fsi_,
    gff = orfgff_
  ) %*>% overlapMap

  f_faa_ <- f_primary@files@aa.file %>>% from_cache
  t_faa_ <- t_primary@files@aa.file %>>% from_cache

  t_orffaa_ <- as_monad(t_primary@files@orffaa.file %>% from_cache)

  t_transorffaa_ <- as_monad(t_primary@files@transorffaa.file %>% from_cache)

  # - align query protein against target genes in the SI
  aa2aa_ <- rmonad::funnel(
    queseq  = f_faa_,
    tarseq  = t_faa_,
    map     = f_si_map_,
    queries = queries
  ) %*>% align_by_map

  # # Run the exact test above, but with the query indices scrambled. For a
  # # reasonably sized genome, there should be very few hits (I should do the
  # # math and make a test)
  # rand_aa2aa_ <- rmonad::funnel(
  #   queseq  = f_faa_,
  #   tarseq  = t_faa_,
  #   map     = f_si_map_,
  #   queries = queries
  # ) %*>% align_by_map(permute=TRUE)

  # - align query protein against ORFs in genomic intervals in SI
  aa2orf_ <- rmonad::funnel(
    queseq  = f_faa_,
    tarseq  = t_orffaa_,
    map     = f_si_map_orf_,
    queries = queries
  ) %*>% align_by_map

  # rand_aa2orf_ <- rmonad::funnel(
  #   queseq  = f_faa_,
  #   tarseq  = t_orffaa_,
  #   map     = f_si_map_orf_,
  #   queries = queries
  # ) %*>% align_by_map(permute=TRUE)

  # - align query protein against ORFs in transcripts in the SI

  # - align query DNA sequence against the SI

  rmonad::funnel(
    si           = si_,
    flag_summary = synder_flags_summary_,
    unassembled  = unassembled_,
    scrambled    = scrambled_,
    indels       = indels_,
    f_si_map     = f_si_map_,
    r_si_map     = r_si_map_,
    aa2aa        = aa2aa_,
    aa2orf       = aa2orf_,
    gapped       = gapped_
  )

}

secondary_data <- function(primary_input, con){

  "
  Loop over all target species, produce the secondary data needed to classify
  orphans. This step will take a long time (many CPU hours).
  "

  focal_species <- con@input@focal_species
  f_primary <- primary_input@species[[focal_species]]

  focal_gff_ <- get_transcripts_from_txdb(f_primary) %>_% {

    "Assert query ids match the names in the GFF file"

    queries <- primary_input@queries
    txnames <- GenomicRanges::mcols(.)$attr
    matches <- queries %in% txnames

    if(!all(matches)){
      msg <- "%s of %s query ids are missing in the GFF. Here is the head: [%s]"
      stop(sprintf(
        msg,
        sum(!matches),
        length(queries),
        paste(head(queries[!matches]), collapse=", ")
      ))
    }
  }

  rmonad::funnel(
    prim = primary_input,
    fgff = focal_gff_
  ) %*>% {

    species <- names(prim@species)
    target_species <- setdiff(species, focal_species)

    ss <- target_species %>%
      lapply(
        function(spec){
          compare_target_to_focal(
            species   = spec,
            fgff      = fgff,
            f_primary = f_primary,
            t_primary = prim@species[[spec]],
            synmap    = prim@synmaps[[spec]],
            queries   = prim@queries,
            con       = con
          )
        }
      )

    names(ss) <- target_species
    rmonad::combine(ss)

  }

}
