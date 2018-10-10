overlapMap <- function(si, gff, type=NULL){

  "Find all features in the gff GenomicRanges object that overlap a search interval."

  if(! class(gff) == "GFF"){
    stop("Expected 'GFF', got '", class(gff), "'")
  }

  # A Hits object
  # from(o) accesses the si indices
  # to(o) acceses the ft indices
  GenomicRanges::findOverlaps(CNEr::second(si), gff) %>% {
    data.frame(
      query            = GenomicRanges::mcols(si[S4Vectors::from(.)])$attr,
      target           = GenomicRanges::mcols(gff[S4Vectors::to(.)])$attr,
      type             = GenomicRanges::mcols(gff[S4Vectors::to(.)])$type,
      qid              = S4Vectors::from(.),
      tid              = S4Vectors::to(.),
      stringsAsFactors = FALSE
    )
  }

}
