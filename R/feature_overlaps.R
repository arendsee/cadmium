overlapMap <- function(si, gff){

  "Find all features in the gff GenomicRanges object that overlap a search interval."

  # A Hits object
  # from(o) accesses the si indices
  # to(o) acceses the ft indices
  GenomicRanges::findOverlaps(CNEr::second(si), gff) %>% {
    data.frame(
      query            = GenomicRanges::mcols(si[from(.)])$attr,
      target           = GenomicRanges::mcols(gff[to(.)])$attr,
      type             = GenomicRanges::mcols(gff[to(.)])$type,
      qid              = from(.),
      tid              = to(.),
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
