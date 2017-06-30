#' Load the configuration file for the Fagin workflow
#'
#' @export
#' @param filename YAML configuration file
load_config <- function(filename="inst/etc/config-template.yml"){
  conf <- yaml::yaml.load_file(filename)
}

#' Thresholds for alignment significance
#'
#' For details on inputs see main package documentation
#'
#' @slot prot2prot      post-correction p-value threshold
#' @slot prot2allorf    post-correction p-value threshold
#' @slot prot2transorf  post-correction p-value threshold
#' @slot dna2dna        post-correction p-value threshold
config_alignment_thresholds <- setClass(
  "config_alignment_thresholds",
  representation(
    prot2prot     = "numeric",
    prot2allorf   = "numeric",
    prot2transorf = "numeric",
    dna2dna       = "numeric"
  ),
  prototype(
    prot2prot     = 0.05,
    prot2allorf   = 0.05,
    prot2transorf = 0.05,
    dna2dna       = 0.05
  )
)

#' Simulation parameters
#'
#' For details on inputs see main package documentation
#'
#' @slot  prot2prot      sample size
#' @slot  prot2allorf    sample size
#' @slot  prot2transorf  sample size
config_alignment_simulation <- setClass(
  "config_alignment_simulation",
  representation(
     prot2prot     = "integer",
     prot2allorf   = "integer",
     prot2transorf = "integer"
  ),
  prototype(
     prot2prot     = 1000L,
     prot2allorf   = 1000L,
     prot2transorf = 1000L
  )
)

#' Alignment settings and thresholds
#'
#' For details on inputs see main package documentation
#'
#' @slot thresholds       p-value thresholds for hit significanc
#' @slot simulation       simulation parameters
#' @slot dna2dna_maxspace largest size of string to compare
#' @slot indel_threshold  I don't remember what this does
config_alignment <- setClass(
  "config_alignment",
  representation(
    thresholds       = "config_alignment_thresholds",
    simulation       = "config_alignment_simulation",
    dna2dna_maxspace = "integer",
    indel_threshold  = "numeric"
  ),
  prototype(
    thresholds = config_alignment_thresholds(),
    simulation = config_alignment_simulation()
  )
)

#' Paths to inputs and the focal species name
#'
#' For details on input requirements, see the main package documentation
#'
#' @slot gff_dir           directory of GFF files (.gff or gff3 extensions)
#' @slot fna_dir           directory of genome files (fasta format)
#' @slot syn_dir           directory of synteny maps
#' @slot tree              directory of phylogenetic tree
#' @slot focal_species     name of the focal species
#' @slot query_gene_list   file containing query genes
config_input <- setClass(
  "config_input",
  representation(
    gff_dir         = "character",
    fna_dir         = "character",
    syn_dir         = "character",
    tree            = "character",
    focal_species   = "character",
    query_gene_list = "character"
  ),
  prototype(
    gff_dir         = "INPUT/gff",
    fna_dir         = "INPUT/fna",
    syn_dir         = "INPUT/syn",
    tree            = "INPUT/tree",
    focal_species   = "Saccharomyces cerevisiae",
    query_gene_list = "INPUT/orphan-list.txt"
  )
)

#' Settings for Synder
#'
#' For details on inputs see main package documentation
#'
#' @slot offsets Base offsets for synder
#' @slot k       How many interlopers shall we allow?
config_synder <- setClass(
  "config_synder",
  representation(
    offsets = "character",
    k       = "integer"
  ),
  prototype(
    offsets = "0111",
    k       = 10L
  )
)

#' The top configuration class
#'
#' The main thing you will need to change is the input.
#'
#' @slot input     Paths to inputs and the focal species name
#' @slot synder    Parameters for synder
#' @slot alignment Alignment configurations
config <- setClass(
  "config",
  representation(
    input     = "config_input",
    synder    = "config_synder",
    alignment = "config_alignment"
  ),
  prototype(
    input     = config_input(),
    synder    = config_synder(),
    alignment = config_alignment()
  )
)
