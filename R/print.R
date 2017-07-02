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


#' @rdname fagin_printer
#' @export 
print.numeric_summary <- function(x, ...){
  cat(sprintf('A numeric vector of length %d\n', x@n))
  print(c(
    "min"      = x@min    ,
                 x@q25    ,
    "median"   = x@median ,
                 x@q75    ,
    "max"      = x@max    ,
    "mean"     = x@mean   ,
    "sd"       = x@sd
  ))
}
setMethod("show", "numeric_summary",
  function(object) plot(object@density)
)

print.synmap_summary <- function(x, ...){
  cat(sprintf("A synmap of %d links\n", x@nrow))
  cat("Width summary:\n")
  print(x@width)
  cat("Score summary:\n")
  print(x@score)
  cat("Query to target log2 ratio\n")
  print(x@query_target_log2_ratio)
}
setMethod("show", "synmap_summary",
  function(object) print(object)
)

# print.seq_summary <- function(x, ...){
#     # seqids           "character",
#     # sizes            "integer",
#     # nseq             "integer",
#     # comp             "matrix",
#     # n_initial_start  "integer",
#     # n_terminal_stop  "integer",
#     # lengths          "integer"
# }
#
# print.dna_summary <- function(x, ...){
#     # n_triple = "integer"
#     # seq
# }
#
# print.faa_summary <- function(x, ...){
#   cat("faa_summary\n")
#     # class_composition  # "numeric"
#   # contains = "seq_summary"
# }
#
# print.gff_summary <- function(x, ...){
#   cat("gff_summary\n")
#     # seqstats     # "data.frame",
#     # mRNA_length  # "numeric_summary",
#     # CDS_length   # "numeric_summary",
#     # exon_length  # "numeric_summary"
# }
#
# print.species_data_files <- function(x, ...){
#   cat("species_data_files\n")
#      # gff.file          # "character",
#      # dna.file          # "character",
#      # aa.file           # "character",
#      # trans.file        # "character",
#      # orfgff.file       # "character",
#      # orffaa.file       # "character",
#      # transorfgff.file  # "character",
#      # transorffaa.file  # "character",
#      # nstring.file      # "character"
# }
#
# print.species_summaries <- function(x, ...){
#   cat("species_summaries\n")
#     # gff.summary          # "gff_summary",
#     # dna.summary          # "dna_summary",
#     # aa.summary           # "faa_summary",
#     # trans.summary        # "dna_summary",
#     # orfgff.summary       # "gff_summary",
#     # orffaa.summary       # "faa_summary",
#     # transorfgff.summary  # "gff_summary",
#     # transorffaa.summary  # "faa_summary",
#     # nstring.summary      # "numeric_summary"
# }
#
# print.species_meta <- function(x, ...){
#   cat("species_meta\n")
#     # files      # "species_data_files",
#     # summaries  # "species_summaries"
# }
#
#
# print.synteny_meta <- function(x, ...){
#   cat("synteny_meta\n")
#     # synmap.file     # "character",
#     # synmap.summary  # "synmap_summary"
# }
#
#
# print.derived_input <- function(x, ...){
#   cat("derived_input\n")
#     # tree           # "phylo",
#     # focal_species  # "character",
#     # queries        # "character",
#     # species        # "list",
#     # synmaps        # "list"
# }
