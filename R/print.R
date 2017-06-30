#' Special printers
#'
#' @param x An object to print
#' @param ... Additional arguments going to God knows where
#' @name fagin_printer
NULL

# local function for printing a field indented with slot depth
prettyCat <- function(tag, value, indent){
  space <- paste0(rep(" ", indent), collapse="")
  cat(sprintf("%s%s = %s\n", space, tag, value))
}

#' @rdname fagin_printer
#' @export 
print.config_alignment_thresholds <- function(x, ...){
  prettyCat("prot2prot",     x@prot2prot,     4)
  prettyCat("prot2allorf",   x@prot2allorf,   4)
  prettyCat("prot2transorf", x@prot2transorf, 4)
  prettyCat("dna2dna",       x@dna2dna,       4)
}

#' @rdname fagin_printer
#' @export 
print.config_alignment_simulation <- function(x, ...){
  prettyCat("prot2prot",     x@prot2prot,     4)
  prettyCat("prot2allorf",   x@prot2allorf,   4)
  prettyCat("prot2transorf", x@prot2transorf, 4)
}

#' @rdname fagin_printer
#' @export 
print.config_alignment <- function(x, ...){
  prettyCat("dna2dna_maxspace", x@dna2dna_maxspace, 2)
  prettyCat("indel_threshold",  x@indel_threshold,  2)
  cat('  Slot "thresholds":\n')
  print(x@thresholds)
  cat('  Slot "simulation":\n')
  print(x@simulation)
}

#' @rdname fagin_printer
#' @export 
print.config_input <- function(x, ...){
  prettyCat("gff_dir"         , x@gff_dir         , 2)
  prettyCat("fna_dir"         , x@fna_dir         , 2)
  prettyCat("syn_dir"         , x@syn_dir         , 2)
  prettyCat("tree"            , x@tree            , 2)
  prettyCat("focal_species"   , x@focal_species   , 2)
  prettyCat("query_gene_list" , x@query_gene_list , 2)
}

#' @rdname fagin_printer
#' @export 
print.config_synder <- function(x, ...){
  prettyCat("offsets" , x@offsets , 2)
  prettyCat("k"       , x@k       , 2)
}

#' @rdname fagin_printer
#' @export 
print.config <- function(x, ...){
  print(x@input)
  print(x@synder)
  print(x@alignment)
}
