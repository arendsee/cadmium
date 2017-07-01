# NOTE: this can fail
extract_tags <- function(attr, tags, get_naked=FALSE, infer_id=FALSE){

  # add ID to tag list if we need to infer ID
  if(infer_id && (! "ID" %in% tags)){
    tags <- c("ID", tags)
  }
  if((get_naked || infer_id) && (! ".untagged" %in% tags)){
    tags <- c(tags, ".untagged")
  }

  pairs <- lapply(stringr::str_split(attr, ";"), stringr::str_split, "=")


  if(length(tags) == 0){
    stop("Error in extract_tags: no tags selected for extraction")
  }

  # annotate untagged
  pairs <- lapply(pairs,
    function(x) {
      lapply(x, function(y) if(length(y) == 1) { c(".untagged", y) } else { y } )
    }
  )

  # This is will appended to the final output table
  ntags <- sapply(pairs, length)

  # assert that each attribute has only a single '='
  lapply(pairs, lapply, length) %>% unlist %>% unlist %>%
    `>`(2) %>% any(na.rm=TRUE) %>%
    'if'(stop("GffError: Illegal second '=' in attribute"))

  # assert that all tags appear only once
  lapply(pairs, lapply, '[', 1) %>%
    lapply(unlist) %>%
    lapply(duplicated) %>% unlist %>% any(na.rm=TRUE) %>%
    'if'(stop("GffError: Repeated tag in 9th column"))

  # assert extracted tags have no commas
  lapply(pairs, lapply, grepl, pattern=",") %>% unlist %>% unlist %>%
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

  # Counter the autocasting horrors of R core
  for(i in seq_len(ncol(a))){
    a[[i]] <- as.character(a[[i]])
  }

  names(a) <- tags
  a$.n_tags <- ntags


  if(infer_id){
    # handle the excrement of AUGUSTUS
    #   * interpret Untagged as ID if no ID is given and if no other tags are present
    a$ID <- ifelse(
      is.na(a$ID)            &  # is this feature has no ID attribute
        (!is.na(a$.untagged)) & # but it does have an attribute with no tag
        a$.n_tags == 1,         # and if this untagged attribute is the only attribute
      a$.untagged,  # if so, assign the untagged attribute to ID
      a$ID          # if not, just use the current ID
    )
  }

  if(!get_naked){
    a$.untagged <- NULL
  }

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
  meta <- extract_tags(gff$attr, tags=c("ID", "Name", "Parent"), infer_id=TRUE)

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
