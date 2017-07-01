# NOTE: this can fail
extract_tags <- function(attr, tags, get_naked=FALSE){
  pairs <- lapply(stringr::str_split(attr, ";"), stringr::str_split, "=")

  if(get_naked) {
    tags <- c(tags, "Untagged")
  }

  # annotate untagged
  pairs <- lapply(pairs,
    function(x) {
      lapply(x, function(y) if(length(y) == 1) { c("Untagged", y) } else { y } )
    }
  )

  # This is will appended to the final output table
  ntags <- sapply(pairs, length)

  # assert that each attribute has only a single '='
  lapply(pairs, lapply, length) %>% unlist %>% unlist %>%
    magrittr::equals(2) %>%
    all(na.rm=TRUE)     %>%
    magrittr::not       %>%
    'if'(stop("GffError: Illegal second '=' in attribute"))

  # assert that all tags appear only once
  lapply(pairs, lapply, '[', 1) %>%
    lapply(unlist) %>%
    sapply(duplicated) %>% any(na.rm=TRUE) %>%
    'if'(stop("GffError: Repeated tag in 9th column"))

  # assert extracted tags have no commas
  sapply(pairs, sapply, grepl, pattern=",") %>%
    any(na.rm=TRUE) %>%
    'if'(stop("GffError: Illegal comma in 9th column"))

  # TODO: Holy Shit!!!!
  gettag <- function(tag) {
    sapply(
      pairs,
      function(attr){
        poop <- sapply(attr, function(x) if(x[1] == tag) { x[2] } else { NA })
        if(all(is.na(poop))){
          NA
        } else {
          poop[!is.na(poop)]
        }
      }
    )
  }

  a <- do.call(cbind, lapply(tags, gettag)) %>% as.data.frame(stringsAsFactors=FALSE)
  names(a) <- tags
  a$.n_tags <- ntags
  a
}

MakeGI <- function(starts, stops, scaffolds, strands=NULL, metadata=NULL, seqinfo=NULL){
  if(is.null(strands)){
    strands=rep("*", length(starts))
  } else {
    strands <- gsub('\\.', '*', strands)
  }
  g <- GenomicRanges::GRanges(
    seqnames = scaffolds,
    ranges   = IRanges::IRanges(starts, stops),
    strand   = strands,
    seqinfo  = seqinfo
  )
  if(!is.null(metadata)){
    GenomicRanges::mcols(g) <- metadata
  }
  g
}


make_GI_with_parent_child_relations <- function(gff){
  meta <- extract_tags(gff$attr, tags=c("ID", "Name", "Parent"), get_naked=TRUE)

  # handle the excrement of AUGUSTUS
  #   * interpret Untagged as ID if no ID is given and if no other tags are present
  meta$ID <- ifelse(
    is.na(meta$ID)            & # is this feature has no ID attribute
      (!is.na(meta$Untagged)) & # but it does have an attribute with no tag
      meta$.n_tags == 1,        # and if this untagged attribute is the only attribute
    meta$Untagged,   # if so, assign the untagged attribute to ID
    meta$ID          # if not, just use the current ID
  )

  meta <- dplyr::select(
    meta,
    id     = .data$ID,
    name   = .data$Name,
    parent = .data$Parent
  )

  meta$type <- gff$type

  # TODO: assert all mRNA and gene entries have an ID
  # TODO: assert all CDS and exon entries have a parent
  # TODO: assert all mRNAs have a parent IFF gene features are present

  # TODO: pass seqinfo, holds sequence lengths and species name
  MakeGI(
    starts    = gff$start,
    stops     = gff$stop,
    scaffolds = gff$seqid,
    strands   = gff$strand,
    metadata  = meta
  ) 
}
