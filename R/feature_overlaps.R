overlapMap <- function(si, gff, type=NULL){

  "Find all features in the gff GenomicRanges object that overlap a search interval."

  stopifnot(class(gff) == 'GFF')

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
