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
#' @slot query_gene_list   file containing query genes (e.g. orphan candidates)
#' @slot control_gene_list file containing control genes
config_input <- setClass(
  "config_input",
  representation(
    gff_dir           = "character",
    fna_dir           = "character",
    syn_dir           = "character",
    tree              = "character",
    focal_species     = "character",
    query_gene_list   = "character",
    control_gene_list = "character"
  ),
  prototype(
    gff_dir           = system.file("yeast", "gff", package='fagin'),
    fna_dir           = system.file("yeast", "fna", package='fagin'),
    syn_dir           = system.file("yeast", "syn", package='fagin'),
    tree              = system.file("yeast", "tree", package='fagin'),
    focal_species     = "Saccharomyces_cerevisiae",
    query_gene_list   = system.file("yeast", "orphan-list.txt", package='fagin'),
    control_gene_list = system.file("yeast", "control-list.txt", package='fagin')
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
    offsets = "integer",
    k       = "integer",
    r       = "numeric",
    trans   = "character"
  ),
  prototype(
    offsets = c(1L,1L),
    k       = 0L,
    r       = 0,
    trans   = "i"
  )
)

#' Settings for Open Reading Frame finding
#'
#' @slot start Start codon (e.g. "ATG|TTG|CTG")
#' @slot stop Stop codon (e.g. "TAA|TGA|TAG")
#' @slot minlen Minimum ORF length (not including start and stop codons)
config_orf <- setClass(
  "config_orf",
  representation(
    start  = "character",
    stop   = "character",
    minlen = "integer"
  ),
  prototype(
    start  = "ATG",
    stop   = "TAA|TGA|TAG",
    minlen = 30L
  )
)

default_decision_tree <- '
gen:
  divisor: Is AA similar to a known gene?
  O1:
    primary: ORFic
    secondary: O1
    label: ORFic - known protein
  trn:
    divisor: Is AA similar to any ORF on known mRNA?
    O2:
      primary: ORFic
      secondary: O2
      label: ORFic - transcribed ORF
    orf:
      divisor: Is AA similar to any ORF anywhere?
      O3:
        primary: ORFic
        secondary: O3
        label: ORFic - unknown ORF
      nuc:
        divisor: Is DNA similar to anything?
        cds:
          divisor: Does match overlap a CDS?
          N1:
            primary: Non-ORFic
            secondary: N1
            label: DNA match to CDS
          exo:
            divisor: Does match overlap an exon?
            N2:
              primary: Non-ORFic
              secondary: N1
              label: SI overlaps exon
            rna:
              divisor: Does match overlap a known transcript?
              N3:
                primary: Non-ORFic
                secondary: N3
                label: DNA match to transcript
              N4:
                primary: Non-ORFic
                secondary: N4
                label: No DNA match to any known gene
        tec:
          U7:
            primary: Unknown
            secondary: U7
            label: skipped
          una:
            divisor: Query maps off target scaffold
            U1:
              primary: Unknown
              secondary: U1
              label: unassembled
            ind:
              divisor: Query maps to zero-length target interval
              U2:
                primary: Unknown
                secondary: U2
                label: possible indel
              nst:
                divisor: Query maps to N-string
                U3:
                  primary: Unknown
                  secondary: U3
                  label: possibly in unknown region
                scr:
                  divisor: Query maps inbetween contiguous block
                  U5:
                    primary: Unknown
                    secondary: U5
                    label: scrambled synteny
                  U6:
                    primary: Unknown
                    secondary: U6
                    label: good syntenic match, no homology
'

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
    input         = "config_input",
    synder        = "config_synder",
    orf           = "config_orf",
    alignment     = "config_alignment",
    decision_tree = "list",
    archive       = "character",
    cache         = "logical"
  ),
  prototype(
    input         = config_input(),
    synder        = config_synder(),
    alignment     = config_alignment(),
    decision_tree = yaml::yaml.load(default_decision_tree),
    archive       = 'ARCHIVE',
    cache         = FALSE
  )
)

#' Get a default configuration object
#'
#' Checks existence of all required files and reports species in the tree
#'
#' @export
config <- function(){
  con <- fagin_config()
  con@input@focal_species <- gsub(" ", "_", con@input@focal_species)
  con
}

#' Validate a configuration
#'
#' @export
validate_config <- function(con){
  if(! file.exists(con@input@tree)){
    stop("Tree file not found")
  }
  if(! file.exists(con@input@gff_dir)){
    stop("GFF directory not found")
  }
  if(! file.exists(con@input@fna_dir)){
    stop("Genome sequence directory not found")
  }
  if(! file.exists(con@input@syn_dir)){
    stop("Synteny map directory not found")
  }
  if(! file.exists(con@input@query_gene_list)){
    stop("Query gene list not found")
  }
  if(! file.exists(con@input@control_gene_list)){
    stop("Control gene list not found")
  }
  # check whether the tree can be loaded
  tree <- ape::read.tree(con@input@tree)
  if(! (con@input@focal_species %in% tree$tip.label)){
    stop("Focal species is not in tree")
  }
  cat("Found tree with species:", paste(tree$tip.label, collapse=", "), "\n")
  cat("Everything looks good so far\n")
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
#' @slot table data.frame of sequence names and lengths
#' @slot comp  matrix base composition
seq_summary <- setClass(
  "seq_summary",
  representation(
    table = "data.frame",
    comp  = "matrix"
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
    n_triple      = "integer",
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
    has_internal_stop = "logical",
    final_residue   = "integer"
  ),
  contains = "seq_summary"
)

#' CDS phase summary 
#'
#' @slot initial_residue   integer
#' @slot final_residue     integer
phase_summary <- setClass(
  "phase_summary",
  representation(
    table = "table",
    incomplete_models = "character"
  )
)

#' Summary of a GFF file
#'
#' @slot table       data.frame
#' @slot mRNA_length numeric_summary
#' @slot CDS_length  numeric_summary
#' @slot exon_length numeric_summary
gff_summary <- setClass(
  "gff_summary",
  representation(
    table       = "data.frame",
    mRNA_length = "numeric_summary",
    CDS_length  = "numeric_summary",
    exon_length = "numeric_summary"
  )
)

#' Summary of a IRanges or GRanges file (without the types of a GFF)
#'
#' @slot table  data.frame
#' @slot width  numeric_summary
granges_summary <- setClass(
  "granges_summary",
  representation(
    table = "data.frame",
    width = "numeric_summary"
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
#' @slot specsum.file     character
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
     nstring.file     = "character",
     specsum.file     = "character"
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
    nstring.summary     = "numeric",
    model_phases        = "list"
  )
)

#' Data summaries and references to full data for a given species
#'
#' @slot files     species_data_files RData files containing the full data
#' @slot summaries character          Filename for saved detailed summaries of all data
#' @slot seqinfo   Seqinfo            Genomic info
species_meta <- setClass(
  "species_meta",
  representation(
    files     = "species_data_files",
    summaries = "character",
    seqinfo   = "Seqinfo"
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
#' @slot queries       A character vector of gene ids (orphan candidates, or whatever)
#' @slot control       A character vector of gene ids (control)
#' @slot species       A list of species_meta objects
#' @slot synteny_maps  A list of synteny maps
derived_input <- setClass(
  "derived_input",
  representation(
    tree          = "phylo",
    focal_species = "character",
    queries       = "character",
    control       = "character",
    species       = "list",
    synmaps       = "list"
  )
)
