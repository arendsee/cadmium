#' Initialize the configuration for the yeast case study
#'
#' @return Config object for yeast using data included in the package 
#' @export
yeast_config <- function(){
  con <- config() 
  con@synder@offsets = c(1L,1L) # offsets for mummer
  con@synder@trans = "p" # percent identity transform mummer 
  con@alignment@dna2dna_maxspace = 1e8L
  con@input@focal_species = "Saccharomyces_cerevisiae"
  con@input@gff <- list(
    "Saccharomyces_arboricola"   = system.file("yeast", "gff", "Saccharomyces_arboricola.gff",   package="fagin")
  , "Saccharomyces_cerevisiae"   = system.file("yeast", "gff", "Saccharomyces_cerevisiae.gff",   package="fagin")
  , "Saccharomyces_eubayanus"    = system.file("yeast", "gff", "Saccharomyces_eubayanus.gff",    package="fagin")
  , "Saccharomyces_kudriavzevii" = system.file("yeast", "gff", "Saccharomyces_kudriavzevii.gff", package="fagin")
  , "Saccharomyces_mikatae"      = system.file("yeast", "gff", "Saccharomyces_mikatae.gff",      package="fagin")
  , "Saccharomyces_paradoxus"    = system.file("yeast", "gff", "Saccharomyces_paradoxus.gff",    package="fagin")
  , "Saccharomyces_uvarum"       = system.file("yeast", "gff", "Saccharomyces_uvarum.gff",       package="fagin")
  )
  con@input@fna <- list(
    "Saccharomyces_arboricola"   = system.file("yeast", "fna", "Saccharomyces_arboricola.fna",   package="fagin")
  , "Saccharomyces_cerevisiae"   = system.file("yeast", "fna", "Saccharomyces_cerevisiae.fna",   package="fagin")
  , "Saccharomyces_eubayanus"    = system.file("yeast", "fna", "Saccharomyces_eubayanus.fna",    package="fagin")
  , "Saccharomyces_kudriavzevii" = system.file("yeast", "fna", "Saccharomyces_kudriavzevii.fna", package="fagin")
  , "Saccharomyces_mikatae"      = system.file("yeast", "fna", "Saccharomyces_mikatae.fna",      package="fagin")
  , "Saccharomyces_paradoxus"    = system.file("yeast", "fna", "Saccharomyces_paradoxus.fna",    package="fagin")
  , "Saccharomyces_uvarum"       = system.file("yeast", "fna", "Saccharomyces_uvarum.fna",       package="fagin")
  )
  con@input@syn <- list(
    "Saccharomyces_arboricola"   = system.file("yeast", "syn", "Saccharomyces_cerevisiae.vs.Saccharomyces_arboricola.syn",   package="fagin")
  , "Saccharomyces_eubayanus"    = system.file("yeast", "syn", "Saccharomyces_cerevisiae.vs.Saccharomyces_eubayanus.syn",    package="fagin")
  , "Saccharomyces_kudriavzevii" = system.file("yeast", "syn", "Saccharomyces_cerevisiae.vs.Saccharomyces_kudriavzevii.syn", package="fagin")
  , "Saccharomyces_mikatae"      = system.file("yeast", "syn", "Saccharomyces_cerevisiae.vs.Saccharomyces_mikatae.syn",      package="fagin")
  , "Saccharomyces_paradoxus"    = system.file("yeast", "syn", "Saccharomyces_cerevisiae.vs.Saccharomyces_paradoxus.syn",    package="fagin")
  , "Saccharomyces_uvarum"       = system.file("yeast", "syn", "Saccharomyces_cerevisiae.vs.Saccharomyces_uvarum.syn",       package="fagin")
  )
  con@input@tree <- system.file("yeast", "tree", package="fagin") 
  con@input@query_gene_list <- system.file("yeast", "orphan-list.txt", package="fagin") 
  con@input@control_gene_list <- system.file("yeast", "control-list.txt", package="fagin") 
  con
}
