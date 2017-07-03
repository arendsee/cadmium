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

# a function for avoiding giant dumps of tabular things
prettyTable <- function(x){
  if(nrow(x) < 40){
    print(x)
  } else {
    print(head(x, 5))
    cat("...\n")
    print(tail(x, 5))
  }
}
prettyVector <- function(x){
  if(length(x) < 400){
    print(x)
  } else {
    print(head(x, 50))
    cat("...\n")
    print(tail(x, 50))
  }
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

print.seq_summary <- function(x, ...){
  cat('Slot "table":\n')
  prettyTable(x@table)
  cat('Slot "comp":\n')
  prettyTable(x@comp)
}
setMethod("show", "seq_summary",
  function(object) print(object)
)

print.dna_summary <- function(x, ...){
  cat('Summary of DNA set with', nrow(x@table) , 'sequence\n')
  cat('Slot "n_triple"\n')
    print(x@n_triple)
  cat('Slot "initial_codon"\n')
    print(x@initial_codon)
  cat('Slot "final_codon"\n')
    print(x@final_codon)
  cat('Slot "table":\n')
    prettyTable(x@table)
  cat('Slot "comp" (%):\n')
    print((100 *colSums(x@comp) / sum(colSums(x@comp))) %>% signif(5))
}
setMethod("show", "dna_summary",
  function(object) print(object)
)

print.faa_summary <- function(x, ...){
  cat('Summary of AA set with', nrow(x@table) , 'sequence\n')
  cat('Slot "initial_residue"\n')
    print(x@initial_residue)
  cat('Slot "initial_residue"\n')
    print(x@final_residue)
  cat('Slot "final_residue"\n')
    prettyTable(x@table)
  cat('Slot "comp"  (%):\n')
    print((100 * colSums(x@comp) / sum(colSums(x@comp))) %>% signif(5))
}
setMethod("show", "faa_summary",
  function(object) print(object)
)

print.gff_summary <- function(x, ...){
  cat('Summary of a "gff_summary" object\n')
  cat('Slot "seqstats"\n')
  prettyTable(x@seqstats)
  cat('Slot "mRNA_length"\n')
  print(x@mRNA_length)
  cat('Slot "CDS_length"\n')
  print(x@CDS_length)
  cat('Slot "exon_length"\n')
  print(x@exon_length)
}
setMethod("show", "gff_summary",
  function(object) print(object)
)

print.species_data_files <- function(x, ...){
  cat('Cached files in "species_data_files" object\n')
  cat(sprintf("gff:         %s\n", x@gff.file))
  cat(sprintf("dna:         %s\n", x@dna.file))
  cat(sprintf("aa:          %s\n", x@aa.file))
  cat(sprintf("trans:       %s\n", x@trans.file))
  cat(sprintf("orfgff:      %s\n", x@orfgff.file))
  cat(sprintf("orffaa:      %s\n", x@orffaa.file))
  cat(sprintf("transorfgff: %s\n", x@transorfgff.file))
  cat(sprintf("transorffaa: %s\n", x@transorffaa.file))
  cat(sprintf("nstring:     %s\n", x@nstring.file))
  cat(sprintf("specsum:     %s\n", x@specsum.file))
}
setMethod("show", "species_data_files",
  function(object) print(object)
)

print.species_summaries <- function(x, ...){
  cat('Collection of data summaries\n')
  cat('\n-------------------------------------------------------------\n')
  cat('Summary of gene model GFF\n\n')
  print(x@gff.summary)
  cat('\n-------------------------------------------------------------\n')
  cat('Summary of genome\n\n')
  print(x@dna.summary)
  cat('\n-------------------------------------------------------------\n')
  cat('Summary of proteins derived from gene model GFF\n\n')
  print(x@aa.summary)
  cat('\n-------------------------------------------------------------\n')
  cat('Summary of mRNA transcripts derived from gene model GFF\n\n')
  print(x@trans.summary)
  cat('\n-------------------------------------------------------------\n')
  cat('Summary of ORFs extracted (naively) from the genome\n\n')
  print(x@orfgff.summary)
  cat('\n-------------------------------------------------------------\n')
  cat('Summary of the translated products of the ORFs\n\n')
  print(x@orffaa.summary)
  cat('\n-------------------------------------------------------------\n')
  cat('Summary of ORFs extracted (naively) from mRNA transcripts\n\n')
  print(x@transorfgff.summary)
  cat('\n-------------------------------------------------------------\n')
  cat('Summary of translated products of the mRNA ORFs\n\n')
  print(x@transorffaa.summary)
  cat('\n-------------------------------------------------------------\n')
  cat('Summary of the strings of Ns in the genome\n\n')
  print(x@nstring.summary)
}
setMethod("show", "species_summaries",
  function(object) print(object)
)

print.species_meta <- function(x, ...){
  cat('Container for all data related to a species\n')
  cat('"@files" contains a list of stored data')
  cat('"@summaries" contains a summary of every record')
}
setMethod("show", "species_meta",
  function(object) print(object)
)

print.derived_input <- function(x, ...){
  cat('Container for all data for all species in this analysis\n')
  cat('Focal species: ', x@focal_species, '\n')
  cat('@tree - phylogenetic tree for all species\n')
  cat('@queries - list of query ids\n')
  cat('@species - list of full data for each species\n')
  cat('@synmaps - list of synteny maps between focal and target species\n')
}
setMethod("show", "derived_input",
  function(object) print(object)
)
