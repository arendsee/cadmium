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

load_gene_models <- function(filename, ...){
  load_gff(filename, tags=c("ID", "Name", "Parent"), infer_id=TRUE, ...)
}

#' Load a GFF file as a GenomicFeatures object
#'
#' @export
#' @param file GFF filename
#' @param tags metadata tags to keep
#' @param get_naked parse untagged attributes as ID if they are they only field
#' @param infer_id if ID is missing, try to create one
#' @param Rmonad wrapped GenomicFeatures object
load_gff <- function(file, tags, get_naked=FALSE, infer_id=FALSE, seqinfo_=NULL){

  "
  Load a GFF. There are many ways this can go wrong. Below is a summary of the
  checks and transforms applied here.

  * check GFF types (character, integer or float) 
  * unify type synonyms
  
  "

  raw_gff_ <-
    {

      "Load the raw GFF file. Raise warnings if columns have incorrect types.
      Allow comments ('#') and use '.' to indicate missing data."

      readr::read_tsv(
        file,
        col_names = c(
          "seqid",
          "source",
          "type",
          "start",
          "stop",
          "score",
          "strand",
          "phase",
          "attr"
        ),
        na        = ".",
        comment   = "#",
        col_types = "ccciidcic"
      )

    } %>_% {

      for(col in c("seqid", "type", "start", "stop")){
        if(any(is.na(.[[col]])))
          stop("GFFError: Column '", col, "' may not have missing values")
      }

  } %>>% {

    "
    Unify all type synonyms. Synonymous sets:

    gene : SO:0000704
    mRNA : messenger_RNA | messenger RNA | SO:0000234
    CDS  : coding_sequence | coding sequence | SO:0000316
    exon : SO:0000147

    The SO:XXXXXXX entries are ontology terms
    "

    gene_synonyms <- 'SO:0000704'
    mRNA_synonyms <- c('messenger_RNA', 'messenger RNA', 'SO:0000234')
    CDS_synonyms  <- c('coding_sequence', 'coding sequence', 'SO:0000316')
    exon_synonyms <- 'SO:0000147'

    .$type <- ifelse(.$type %in% gene_synonyms, 'gene', .$type)
    .$type <- ifelse(.$type %in% mRNA_synonyms, 'mRNA', .$type)
    .$type <- ifelse(.$type %in% CDS_synonyms,  'CDS',  .$type)
    .$type <- ifelse(.$type %in% exon_synonyms, 'exon', .$type)

    .

  } %>>% {

    "
    Replace transcript and coding_exon (and their synonyms) with mRNA and exon,
    respectively. This is not formally correct, but is probably the right thing
    to do. Since this is questionable, a warning is emitted if any replacements
    are made.

    The following nearly synonymous sets are merged:

    mRNA : transcript | SO:0000673
    exon : SO:0000147 | coding_exon | coding exon | SO:0000195
    "

    mRNA_near_synonyms <- c('transcript', 'SO:0000673')
    exon_near_synonyms <- c('SO:0000147', 'coding_exon', 'coding exon', 'SO:0000195')

    if(any(.$type %in% mRNA_near_synonyms)){
        .$type <- ifelse(.$type %in% mRNA_near_synonyms, 'mRNA', .$type)
        warning("Substituting transcript types for mRNA types, this is probably OK")
    }

    if(any(.$type %in% exon_near_synonyms)){
        .$type <- ifelse(.$type %in% exon_near_synonyms, 'exon', .$type)
        warning("Substituting transcript types for exon types, this is probably OK")
    }

    .

  }


  tags_ <- tags %v>%
  {

    "
    Internal. Setup and check tag list.

    1) add ID to tag list if we need to infer ID
    2) sets a temporary tag for untagged entry
    3) assert at least one tag is pressent (otherwise nothing would be done)
    "

    if(infer_id && (! "ID" %in% .)){
      . <- c("ID", .)
    }
    if((get_naked || infer_id) && (! ".U" %in% .)){
      . <- c(., ".U")
    }

    if(length(.) == 0){
      stop("No tags selected for extraction")
    }

    .

  }

  raw_gff_ %>>% {

    "Extract the attribute column"

    .[[9]]

  } %>>% {

    "Split attribute column into individual fields; expressed as a dataframe
    with columns [order | ntags | tag | value]. `order` records the original
    ordering of the GFF file (which will be lost). `ntags` is a count of the
    total number of tags for a GFF row; it is a temporary column. Untagged
    values are given the temporary tag '.U', e.g. 'gene01' -> '.U=gene01'."

    tibble::data_frame(
      attr  = stringr::str_split(., ";"),
      order = seq_len(length(.))
    )                                                                                  %>>%
      dplyr::mutate(ntags = sapply(attr, length))                                      %>>%
      tidyr::unnest(attr)                                                              %>>%
      dplyr::mutate(attr = ifelse(grepl('=', attr), attr, paste(".U", attr, sep="="))) %>>%
      tidyr::separate_(
        col   = "attr",
        into  = c("tag", "value"),
        sep   = "=",
        extra = "merge"
      )

  } %>>% rmonad::funnel(tags=tags_) %*>% {

    "Ignore any tags other than the specified ones"

    dplyr::filter(., tag %in% tags)

  } %>_% {

    "Assert there are no commas in the extracted attribute values. These are
    legal according to the GFF spec, but I do not yet support them."

    if(any(grepl(",", .$value))){
      stop("GFFError: commas not supported in attribute tags")
    }

  } %>% rmonad::funnel(tags=tags_) %*>% {

    "Give each tag its own column"

    if(nrow(.) > 0){
      . <- tidyr::spread(., key="tag", value="value")
      for(tag in tags){
        if(!tag %in% names(.)){
          .[[tag]] <- NA_character_
        }
      }
      .
    } else {
      .$tag   = NULL
      .$value = NULL
      for(tag in tags){
        .[[tag]] <- character(0)
      }
      .
    }

  } %>>% {

    "Consider features with the parent '-' to be missing"

    if("Parent" %in% names(.)){
      .$Parent <- ifelse(.$Parent == "-", NA, .$Parent)
    }
    .

  } %>% rmonad::funnel(infer_id=infer_id) %*>% {

    "If no ID is given, but there is one untagged field, and if there are no
    other fields, then cast the untagged field as an ID. This is needed to
    accommodate the reprehensible output of AUGUSTUS."

    if(infer_id && ".U" %in% names(.)){
      # handle the excrement of AUGUSTUS
      #   * interpret .U as ID if no ID is given and if no other tags are present
      .$ID <- ifelse(
        is.na(.$ID)       & # is this feature has no ID attribute
          (!is.na(.$.U))  & # but it does have an attribute with no tag
          .$ntags == 1,     # and if this untagged attribute is the only attribute
        .$.U,         # if so, assign the untagged attribute to ID
        .$ID          # if not, just use the current ID
      )
    }

    .

  } %>% rmonad::funnel(gff=raw_gff_) %*>% {

    "Merge the attribute columns back into the GFF, remove temporary columns."

    gff$order <- seq_len(nrow(gff))

    merge(., gff, all=TRUE) %>%
      dplyr::arrange(order) %>%
      dplyr::select(-.U, -order, -ntags, -attr)

  } %>_% {

    "
    Assert the parent/child relations are correct
    "

    if(all(c("ID", "Parent") %in% names(.))){
      parents <- subset(., type %in% c("CDS", "exon"))$Parent
      parent_types <- subset(., ID %in% parents)$type

      if(any(parent_types == "gene"))
        warning("Found CDS or exon directly inheriting from a gene, this may be fine.")

      if(! all(parent_types %in% c("gene", "mRNA"))){
        offenders <- parent_types[!(parent_types %in% c("gene", "mRNA"))]
        msg <- "Found CDS or exon with unexpected parent: [%s]"
        warning(sprintf(msg, paste0(unique(offenders), collapse=", ")))
      }

      if( any(is.na(parents)) )
        stop("Found CDS or exon with no parent")

      if(! any(duplicated(.$ID, incomparables=NA)))
        warning("IDs are not unique, this is probably bad")
    }

  } %>>% {

    "Load GFF into a GenomicRanges object"

    gi <- MakeGI(
      starts    = .$start,
      stops     = .$stop,
      scaffolds = .$seqid,
      strands   = .$strand,
      metadata  = .[,c('ID','Name','Parent')]
    )

    GenomicRanges::mcols(gi)$type <- .$type

    GenomicRanges::seqinfo(gi) <- seqinfo_

    gi

  } %>>% {

    "
    From the GenomicRanges object, create a transcript database as a TxDb object

    makeTxDb has two required inputs (both as data frames). These notes are from (?makeTxDb):
      1) transcripts - one row per transcript with the following columns
         * tx_id - Transcript ID as an integer vector without NA or duplicates
         * tx_name - [optional] Transcript names as character vector
         * tx_chrm - transcript chromosomes as character vector with no NA
         * tx_strand - transcript strand as factor or character vector with no NA from [+-]
         * tx_start - transcript start position
         * tx_end - transcript end position
      2) splicings - must have one row for each exon in each transcript, each
         row describes the exon and the CDS contained in it (if any) 
         * tx_id
         * exon_rank
         * exon_id
         * exon_name
         * exon_chrom
         * exon_start
         * exon_end
         * cds_id
         * cds_name
         * cds_start
         * cds_end
    makeTxDb also has the following optional inputs
      3) genes - must have one row for every transcript in every gene (often equals 1)
         * tx_id/tx_name - genes must have either tx_id OR tx_name, but not both
         * gene_id - Gene ID, character vector or factor with no NAs
      4) chrominfo - must have 1 row per chromosome
         * chrom - Chromosome name as character vector or factor with no NAs or duplicates
         * length - chromosome length as integer vector with either all or no NAs
         * is_circular - logical vector possibly with NAs
    "

    meta <- GenomicRanges::mcols(.)

    # NOTE: This cannot just be `is_trans <- meta$type == "mRNA"` because some
    # exons are recorded as direct children of a "gene" feature. So I list as a
    # transcript anything that is the parent of an exon or CDS.
    is_trans <- meta$ID %in% meta$Parent[meta$type %in% c("exon", "CDS")]

    GenomicRanges::mcols(.)$type = ifelse(is_trans, "mRNA", meta$type)

    GenomicFeatures::makeTxDbFromGRanges(.)

  }

}
