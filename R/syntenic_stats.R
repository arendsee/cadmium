summarize_syntenic_flags <- function(si, queries){
  # require(reshape2)
  # require(dplyr)

  if(!all(queries %in% si$attr)){
    stop(
      "expected the 'attr' column of the synder search interval table to ",
      "include all ids in the query list. It doesn't. This probably means ",
      "'fagin' and 'synder' are out of sync."
    )
  }

  reliable <- si[ ! si$inbetween, ]$attr %>% unique

  reliable

  # scrambled <- setdiff(si$query$seqid, reliable)
  #
  # d <- data.frame(
  #   seqid     = si$query$seqid,
  #   lo_flag   = si$target$lo_flag %>% as.character %>% as.integer,
  #   hi_flag   = si$target$hi_flag %>% as.character %>% as.integer,
  #   inbetween = si$target$inbetween,
  #   is_orphan = si$query$seqid %in% query$orphan
  # ) %>%
  #   dplyr::group_by(seqid) %>%
  #   dplyr::summarize(
  #     inbetween    = all(inbetween),
  #     lo_bound     = any(lo_flag < 2),
  #     hi_bound     = any(hi_flag < 2),
  #     doubly_bound = any(lo_flag < 2 & hi_flag < 2),
  #     unbound      = all(lo_flag > 1 & hi_flag > 1),
  #     beyond       = any(lo_flag == 4 | hi_flag == 4),
  #     is_orphan    = all(is_orphan)
  #   ) %>%
  #   dplyr::group_by(is_orphan) %>%
  #   dplyr::summarize(
  #     inbetween    = mean(inbetween),
  #     lo_bound     = mean(lo_bound),
  #     hi_bound     = mean(hi_bound),
  #     doubly_bound = mean(doubly_bound),
  #     unbound      = mean(unbound),
  #     beyond       = mean(beyond)
  #   ) %>%
  #   dplyr::select(-matches('is_orphan')) %>%
  #   t %>%
  #   set_colnames(c("not_orphan", "orphan"))
  #
  # list(
  #   scrambled = scrambled,
  #   sum = d
  # )
}
