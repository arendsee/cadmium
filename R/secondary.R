compare_target_to_focal <- function(
  species,
  gff,
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

  synmap_ <- as_monad( from_cache(synmap@synmap.file) )

  si_ <- synmap_ %>>% {
    rmonad::funnel(
      syn =
        rmonad::funnel(
          first=GenomicRanges::GRanges(
            seqnames = .$qseqid,
            ranges   = IRanges(.$qstart, .$qstop),
            strand   = '+',
            seqinfo  = f_primary@seqinfo
          ),
          second=GenomicRanges::GRanges(
            seqnames = .$tseqid,
            ranges   = IRanges(.$tstart, .$tstop),
            strand   = .$strand,
            seqinfo  = t_primary@seqinfo,
            score    = .$score
          )
        ) %*>% CNEr::GRangePairs,
      gff = gff
    )
  } %*>%
  synder::search(
    trans   = con@synder@trans,
    k       = con@synder@k,
    r       = con@synder@r,
    offsets = con@synder@offsets
  )

  # synder_flags_summary_ <- si_ %>>% summarize_synder_flags,
  #
  # unassembled_ <- si_ %>>% find_unassembled,
  #
  # indels_ <- funnel(si_, t_primary) %*>% find_indels,
  #
  # gapped_ <- funnel(si_, t_primary@nstrings) %*>% find_gapped,

  # - find target CDS and mRNA that are in SI for each query

  # - align query protein against target genes in the SI

  # - align query protein against ORFs in transcripts in the SI

  # - align query protein against ORFs in genomic intervals in SI

  # - align query DNA sequence against the SI

  si_

}

secondary_data <- function(primary_input, con){

  "
  Loop over all target species, produce the secondary data needed to classify
  orphans. This step will take a long time (many CPU hours).
  "

  gff_ <- {
    primary_input@species[[con@input@focal_species]] %>%
      {from_cache(.@files@gff.file, type='sqlite')}
  } %>_% {

    "Assert query ids match the names in the GFF file"

    queries <- primary_input@queries
    txnames <- GenomicRanges::mcols(GenomicFeatures::transcripts(.))$tx_name
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

  } %>>%
    GenomicFeatures::transcripts %>>%
    synder::as_gff(id_tag="tx_name")

  rmonad::funnel(
    prim = primary_input,
    gff = gff_
  ) %*>% {

    focal_species <- con@input@focal_species
    species <- names(prim@species)
    target_species <- setdiff(species, focal_species)

    ss <- target_species %>%
      lapply(
        function(spec){
          compare_target_to_focal(
            species   = spec,
            gff       = gff,
            f_primary = prim@species[[focal_species]],
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
