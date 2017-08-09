summarize_syntenic_flags <- function(si, queries){

  met <- GenomicRanges::mcols(si)

  if(!all(queries %in% met$attr)){
    stop(
      "expected the 'attr' column of the synder search interval table to ",
      "include all ids in the query list. It doesn't. This probably means ",
      "'fagin' and 'synder' are out of sync."
    )
  }

  data.frame(
    seqid     = met$attr,
    lo_flag   = met$l_flag,
    hi_flag   = met$r_flag,
    inbetween = met$inbetween,
    is_orphan = met$attr %in% queries
  ) %>%
    dplyr::group_by(seqid) %>%
    dplyr::summarize(
      inbetween    = all(inbetween),
      lo_bound     = any(lo_flag < 2),
      hi_bound     = any(hi_flag < 2),
      doubly_bound = any(lo_flag < 2 & hi_flag < 2),
      unbound      = all(lo_flag > 1 & hi_flag > 1),
      beyond       = any(lo_flag == 4 | hi_flag == 4),
      is_orphan    = all(is_orphan)
    ) %>%
    dplyr::group_by(is_orphan) %>%
    dplyr::summarize(
      inbetween    = mean(inbetween),
      lo_bound     = mean(lo_bound),
      hi_bound     = mean(hi_bound),
      doubly_bound = mean(doubly_bound),
      unbound      = mean(unbound),
      beyond       = mean(beyond),
      total        = n()
    ) %>%
    dplyr::select(-dplyr::matches('is_orphan')) %>%
    t %>%
    set_colnames(c("not_orphan", "orphan"))

}

find_scrambled <- function(si){


  "Get the ids of all queries that map to scrambled regions of the target
  genome."

  met <- GenomicRanges::mcols(si)

  met[ met$l_flag > 1 & met$r_flag > 1 & met$inbetween, ]$attr %>% unique

}

find_unassembled <- function(si){

  "Find all genes that map to regions at the extrema of scaffolds"

  met <- GenomicRanges::mcols(si)

  met[ met$l_flag == 3 | met$r_flag == 3, ]$attr %>% unique

}
