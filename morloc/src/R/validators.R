#' Data types and validators
#'
#' Simple types
#' \itemize{
#'   \item Bool
#'   \item Character
#'   \item Numeric
#'   \item Integer
#'   \item Count
#' }
#'
#' Compound types
#' \itemize{
#'   \item Synmap
#'   \item SearchInterval
#'   \item GFF
#'   \item Genome
#'   \item CDS
#'   \item Protein
#'   \item Phylogeny
#'   \item GeneID
#'   \item ProtMSA
#'   \item NuclMSA
#' }
#'
#' @param x data to be checked
#' @name validator



# --- inputs ----------------------------------------------

#' @rdname validator
#' @export
validate_synmap <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_synder <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_gff <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_fna <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_frn <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_faa <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_species_tree <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_query_list <- function(x){
  TRUE
}


# --- others ----------------------------------------------

#' @rdname validator
#' @export
validate_aa_aln_result <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_trees_labels_summary <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_origins <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_feature_matrix <- function(x){
  TRUE
}


# --- tuples ----------------------------------------------

#' @rdname validator
#' @export
validate_synmap_synder <- function(x){
  TRUE
}

#' @rdname validator
#' @export
validate_gff_fna <- function(x){
  TRUE
}
