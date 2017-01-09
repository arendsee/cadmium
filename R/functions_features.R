#' Feature vector generating functions
#'
#' Is AA similar to a known gene
#' gen
#' primary: ORFic
#' secondary: O1
#' label: ORFic - known protein
#'
#' Is AA similar to any ORF on known mRNA?
#' trn
#' divisor: Is AA similar to any ORF on known mRNA?
#' primary: ORFic
#' secondary: O2
#' label: ORFic - transcribed ORF
#'
#' Is AA similar to any ORF anywhere?
#' orf
#' divisor: Is AA similar to any ORF anywhere?
#' primary: ORFic
#' secondary: O3
#' label: ORFic - unknown ORF
#'
#' Is DNA similar to anything?
#' nuc
#'
#' Does search interval overlap a CDS?
#' cds
#' divisor: Does search interval overlap a CDS?
#' N1:
#' primary: Non-ORFic
#' secondary: N1
#' label: SI overlaps CDS
#'
#' Does search interval overlap an exon?
#' rna
#' divisor: Does search interval overlap an exon?
#' N2:
#' primary: Non-ORFic
#' secondary: N2
#' label: SI overlaps exon
#
#' Was the entry skipped for technical reasons?
#' tec:
#' U7:
#' primary: Unknown
#' secondary: U7
#' label: skipped
#'
#' Query maps off target scaffold
#' una:
#' divisor: Query maps off target scaffold
#' U1:
#' primary: Unknown
#' secondary: U1
#' label: unassembled
#'
#' Query maps to zero-length target interval
#' ind:
#' divisor: Query maps to zero-length target interval
#' U2:
#' primary: Unknown
#' secondary: U2
#' label: possible indel
#'
#' Query maps to N-string
#' nst:
#' divisor: Query maps to N-string
#' U3:
#' primary: Unknown
#' secondary: U3
#' label: possibly in unknown region
#'
#' Query maps to target interval smaller than self
#' res:
#' divisor: Query maps to target interval smaller than self
#' U4:
#' primary: Unknown
#' secondary: U4
#' label: possibly resized
#'
#' Query maps inbetween contiguous block
#' scr:
#' divisor: Query maps inbetween contiguous block
#' U5:
#' primary: Unknown
#' secondary: U5
#' label: scrambled synteny
#'
#' Good syntenic match, no homology
#' U6:
#' primary: Unknown
#' secondary: U6
#'
#' @param x input data
#' @return logical vector
#' @name feature_builders


#' @rdname feature_builders
#' @export
find_gen <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_trn <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_orf <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_nuc <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_cds <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_rna <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_tec <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_una <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_ind <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_nst <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_res <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_scr <- function(x){
  NULL
}

#' @rdname feature_builders
#' @export
find_unknown_unknowns <- function(x){
  NULL
}
