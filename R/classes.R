# ================================================================= #
#                  C O N F I G   C L A S S E S                      #
# ================================================================= #


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
    simulation = config_alignment_simulation(),
    dna2dna_maxspace = 1e7L,
    indel_threshold  = 0.25
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
fagin_config <- setClass(
  "fagin_config",
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

#' Get a default configuration object
#'
#' @export
config <- function(){
  fagin_config()
}


# ================================================================= #
#                    D A T A   C L A S S E S                        #
# ================================================================= #

setOldClass("density")
setOldClass("htest")
setOldClass("phylo")


#' Summary of a numeric vector
#'
#' @slot min              numeric  minimum value
#' @slot q25              numeric  25th quantile
#' @slot median           numeric  median
#' @slot q75              numeric  75th quantile
#' @slot max              numeric  maximum value
#' @slot mean             numeric  mean
#' @slot sd               numeric  standard deviation
#' @slot n                integer  total number of elements
#' @slot density          density  kernel density
numeric_summary <- setClass(
  "numeric_summary",
  representation(
    min     = "numeric",
    q25     = "numeric",
    median  = "numeric",
    q75     = "numeric",
    max     = "numeric",
    mean    = "numeric",
    sd      = "numeric",
    n       = "integer",
    density = "density"
  )
)

#' Summary of a single synteny map
#'
#' @slot nrow   integer
#' @slot width  numeric_summary
#' @slot score  numeric_summary
#' @slot query_target_log2_ratio numeric_summary
synmap_summary <- setClass(
  "synmap_summary",
  representation(
    # length_to_score_coref   = "matrix"
    nrow  = "integer",
    width = "numeric_summary",
    score = "numeric_summary",
    query_target_log2_ratio = "numeric_summary"
  )
)

#' Summary of a sequence
#'
#' This is a parent to both DNA and protien summary classes
#'
#' @slot seqids          character the name of each entry
#' @slot sizes           integer lengths of each sequence
#' @slot nseq            integer total number of sequences
#' @slot comp            matrix base composition
#' @slot comp_dist       list 
#' @slot n_non_canonical integer
#' @slot length_sum      numeric_summary
seq_summary <- setClass(
  "seq_summary",
  representation(
    seqids          = "character",
    sizes           = "integer",
    nseq            = "integer",
    comp            = "matrix",
    lengths         = "integer"
  )
)

#' Summary of a set of DNA sequences
#'
#' @slot n_triple        integer
#' @slot initial_codon  integer
#' @slot final_codon    integer
dna_summary <- setClass(
  "dna_summary",
  representation(
    n_triple       = "integer",
    initial_codon = "integer",
    final_codon   = "integer"
  ),
  contains = "seq_summary"
)

#' Summary of a set of protein sequences
#'
#' @slot initial_residue   integer
#' @slot final_residue     integer
faa_summary <- setClass(
  "faa_summary",
  representation(
    initial_residue = "integer",
    final_residue   = "integer"
  ),
  contains = "seq_summary"
)

#' Summary of a GFF file
#'
#' @slot seqstats    data.frame
#' @slot mRNA_length numeric_summary
#' @slot CDS_length  numeric_summary
#' @slot exon_length numeric_summary
gff_summary <- setClass(
  "gff_summary",
  representation(
    seqstats    = "data.frame",
    mRNA_length = "numeric_summary",
    CDS_length  = "numeric_summary",
    exon_length = "numeric_summary"
  )
)

#' Summary of a IRanges or GRanges file (without the types of a GFF)
#'
#' @slot seqstats  data.frame
#' @slot width     numeric_summary
granges_summary <- setClass(
  "granges_summary",
  representation(
    seqstats = "data.frame",
    width    = "numeric_summary"
  )
)

#' References to each of the RData files for a given species
#'
#' @slot gff.file         character
#' @slot dna.file         character
#' @slot aa.file          character
#' @slot trans.file       character
#' @slot orfgff.file      character
#' @slot orffaa.file      character
#' @slot transorfgff.file character
#' @slot transorffaa.file character
#' @slot nstring.file     character
species_data_files <- setClass(
  "species_data_files",
  representation(
     gff.file         = "character",
     dna.file         = "character",
     aa.file          = "character",
     trans.file       = "character",
     orfgff.file      = "character",
     orffaa.file      = "character",
     transorfgff.file = "character",
     transorffaa.file = "character",
     nstring.file     = "character"
   )
)

#' Summaries of the data stored for each species
#'
#' @slot gff.summary      gff_summary
#' @slot dna.summary      dna_summary
#' @slot aa.summary       faa_summary
#' @slot trans.summary    dna_summary
#' @slot orfgff.summary   gff_summary
#' @slot orffaa.summary   faa_summary
#' @slot transorf.summary dna_summary
#' @slot nstring.summary  numeric_summary
species_summaries <- setClass(
  "species_summaries",
  representation(
    gff.summary         = "gff_summary",
    dna.summary         = "dna_summary",
    aa.summary          = "faa_summary",
    trans.summary       = "dna_summary",
    orfgff.summary      = "granges_summary",
    orffaa.summary      = "faa_summary",
    transorfgff.summary = "granges_summary",
    transorffaa.summary = "faa_summary",
    nstring.summary     = "numeric"
  )
)

#' Data summaries and references to full data for a given species
#'
#' @slot files     species_data_files RData files containing the full data
#' @slot summaries species_summaries  Detailed summaries of all data
species_meta <- setClass(
  "species_meta",
  representation(
    files     = "species_data_files",
    summaries = "species_summaries"
  )
)

#' Summaries of every synteny map and references to full data
#'
#' @slot synmap.file    filename,
#' @slot synmap.summary synmap_summary
synteny_meta <- setClass(
  "synteny_meta",
  representation(
    synmap.file    = "character",
    synmap.summary = "synmap_summary"
  )
)

#' Class containing all data required for a Fagin run
#'
#' @slot tree          Phologenetic tree of all species in the analysis
#' @slot focal_species The focal species
#' @slot queries       A character vector of species ids
#' @slot species       A list of species_meta objects
#' @slot synteny_maps  A list of synteny maps
derived_input <- setClass(
  "derived_input",
  representation(
    tree          = "phylo",
    focal_species = "character",
    queries       = "character",
    species       = "list",
    synmaps       = "list"
  )
)
