overlapMap <- function(si, tgff){

  "Find all features in the tgff GenomicRanges object that overlap a search interval."

  # A Hits object
  # from(o) accesses the si indices
  # to(o) acceses the ft indices
  GenomicRanges::findOverlaps(CNEr::second(si), tgff) %>% {
    data.frame(
      query            = GenomicRanges::mcols(si[from(.)])$attr,
      target           = GenomicRanges::mcols(tgff[to(.)])$attr,
      type             = GenomicRanges::mcols(tgff[to(.)])$type,
      stringsAsFactors = FALSE
    )
  }

}

findQueryGaps <- function(si, nstring){

  # data.frame(
  #   query = sique$seqid[from(over)],
  #   length = nstring[to(over)] %>% width,
  #   stringsAsFactors=FALSE
  # )

}
