# attr <- read.table("z", stringsAsFactors=FALSE)[[1]]

# NOTE: this can fail
extract_tags <- function(attr, tags, get_naked=FALSE, infer_id=FALSE){

  # add ID to tag list if we need to infer ID
  if(infer_id && (! "ID" %in% tags)){
    tags <- c("ID", tags)
  }
  if((get_naked || infer_id) && (! ".U" %in% tags)){
    tags <- c(tags, ".U")
  }

  if(length(tags) == 0){
    stop("Error in extract_tags: no tags selected for extraction")
  }

  # temporary ids,just needed for grouping upon
  input_order = 1:length(attr)

  d <- stringr::str_split(attr, ";")

  # This is will appended to the final output table
  ntags <- sapply(d, length)

  d <- data.frame(
    id = rep.int(input_order, times=ntags),
    # fields with no tag are given the default tag ".U", e.g. "foo" -> ".U=foo"
    fields = unlist(d) %>% {ifelse(grepl("=", .), ., paste0(".U=", .))} ,
    stringsAsFactors=FALSE
  ) %>%
    tidyr::separate_(
      col   = "fields",
      into  = c("tag", "value"),
      sep   = "=",
      extra = "merge"
    ) %>%
    dplyr::filter(.data$tag %in% tags)

  if(any(grepl(",", d$value))){
    stop("GFFError: commas not supported in attribute tags")
  }

  if(nrow(d) > 0){
    d <- tidyr::spread(d, key="tag", value="value")
  } else {
    d$tag = NULL
    d$value = NULL
  }

  # fills in empty rows for empty for features with no remaining tags
  d <- merge(d, data.frame(id=input_order), all=TRUE)

  # All tags should be present in the output table, even those that are not
  # represented in the attribute column. Here I set such tags to columns of NA.
  for(missing_tag in setdiff(tags, names(d))){
    d[[missing_tag]] = NA_character_
  }

  d$.n_tags <- ntags

  if(infer_id && ".U" %in% names(d)){
    # handle the excrement of AUGUSTUS
    #   * interpret .U as ID if no ID is given and if no other tags are present
    d$ID <- ifelse(
      is.na(d$ID)       & # is this feature has no ID attribute
        (!is.na(d$.U))  & # but it does have an attribute with no tag
        d$.n_tags == 1,   # and if this untagged attribute is the only attribute
      d$.U,         # if so, assign the untagged attribute to ID
      d$ID          # if not, just use the current ID
    )
  }

  if(!get_naked){
    d$.U <- NULL
  }
  d$id <- NULL

  d
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
