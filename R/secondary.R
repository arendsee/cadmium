compare_target_to_focal <- function(
  species,
  f_primary,
  t_primary,
  synmap,
  queries,
  config
){

  "
  Collect secondary data for one species. The main output is a feature table.
  Internal data is summarized (for use in downstream diagnostics or
  side-analyses) and cached.
  "

  # # - run synder
  # synder_ <-
  #
  #   rmonad::funnel(
  #     syn= synmap
  #     gff=
  #     tcl=
  #     qcl=
  #     trans=
  #     k=
  #     r=
  #     offsets=
  #   ) %*>%
  #   synder::search

  # - summarize synder flags

  # - identify scrambled intervals

  # - identify unassembled intervals

  # - find indels

  # - map intervals to N-string gaps

  # - find target CDS and mRNA that are in SI for each query

  # - align query protein against target genes in the SI

  # - align query protein against ORFs in transcripts in the SI

  # - align query protein against ORFs in genomic intervals in SI

  # - align query DNA sequence against the SI

  as_monad(species)

}

secondary_data <- function(primary_input, config){

  "
  Loop over all target species, produce the secondary data needed to classify
  orphans. This step will take a long time (many CPU hours).
  "

  primary_input@species %>>% {

    focal_species <- config@input@focal_species

    ss <- setdiff(., focal_species) %>% lapply(
      function(spec){
        compare_target_to_focal(
          species   = spec,
          f_primary = primary_input@species[[focal_species]],
          t_primary = primary_input@species[[spec]],
          synmap    = primary_input@synmaps[[spec]],
          queries   = primary_input@queries,
          config    = config
        )
      }
    )

    names(ss) <- .
    rmonad::combine(ss)

  }

}
