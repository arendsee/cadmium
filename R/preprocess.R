#' Load a nexus format tree
#'
#' @param treefile A filename for a nexus tree
#' @return a phylo object
#' @export
load_tree_io <- function(treefile){
  ape::read.tree(treefile)
}

#' Get species from tree
#'
#' @param tree phylo object representing all under inspection
#' @return character vector of unique species
#' @export
get_species_from_tree <- function(tree){
  tree$tip.label
}

read_queries <- function(filename){
  readr::read_table(filename, col_names=FALSE, col_types="c", comment="#")
}


### Full TODO:
# -- required for release --------------------------------------
# 1. implement the functions below (completes preprocessing)
# 2. add in error handling using monadR
# -- fairly mindless refactoring ---
# 3. classify all tertiary results from the classify.R script
# 4. tie all those pieces together
# -- need to do eventually -------------------------------------
# 5. write test suite
# 6. write report generator
# 7. write vignette

### TODO: fill in the following functions ###########################

to_cache    <- function(x, species_name, tag) { }

# *******************************************************************


load_species <- function(species_name, input){
  # Primary data - required inputs to Fagin
  dna <- load_dna( get_genome_filename ( species_name, dir ))
  gff <- load_gff( get_gff_filename    ( species_name, dir ))

  # Gene model independent data, derived only from genome
  nstring <- derive_nstring(dna)
  orfgff  <- derive_orfgff(dna)
  orffaa  <- derive_orffaa(dna, orfgff)

  # Gene model derived data
  aa          <- derive_aa(dna, gff)
  trans       <- derive_trans(dna, gff)
  transorfgff <- derive_transorfgff(trans)
  transorffaa <- derive_transorffaa(trans)

  specsum <- new("species_summaries",
    gff.summary         = summarize_gff(gff),
    dna.summary         = summarize_fna(dna),
    aa.summary          = summarize_faa(aa),
    trans.summary       = summarize_fna(trans),
    orfgff.summary      = summarize_gff(orfgff),
    orffaa.summary      = summarize_faa(orffaa),
    transorfgff.summary = summarize_gff(transorfgff),
    transorffaa.summary = summarize_faa(transorffaa),
    nstring.summary     = summarize_nstring(nstring)
  )

  specfile <- new("species_data_file",
    gff.file         = to_cache( gff         , species_name , "gff"         ),
    dna.file         = to_cache( dna         , species_name , "dna"         ),
    aa.file          = to_cache( aa          , species_name , "aa"          ),
    trans.file       = to_cache( trans       , species_name , "trans"       ),
    orfgff.file      = to_cache( orfgff      , species_name , "orfgff"      ),
    orffaa.file      = to_cache( orffaa      , species_name , "orffaa"      ),
    transorfgff.file = to_cache( transorfgff , species_name , "transorfgff" ),
    transorffaa.file = to_cache( transorffaa , species_name , "transorffaa" ),
    nstring.file     = to_cache( nstring     , species_name , "nstring"     )
  )

  new("species_meta", files = specfile, summaries = specsum)

}

load_synmaps <- function(target_species, focal_species, syndir){

  # TODO: replace this hard-coded pattern
  filename <- get_synmap_filename(focal_species, target_species)

  synfile <- load_synmap(filename)

  synsum <- summarize_syn(synfile)

  new(
    "synteny_meta",
    synmap.file    = synfile,
    synmap.summary = synsum
  )
}

#' Derive secondary data from the minimal required inputs
#'
#' The required inputs to Fagin are genomes fasta files, GFFs, synteny maps, a
#' species tree, and the focal species name. This function derives all required
#' secondary data, but does no analysis specific to the query intervals (e.g.
#' the orphan genes candidates).
#'
#' This function does not handle the full data validation. All individual
#' pieces are validated, but not the relationships between them. This is the
#' role of the \code{validate_derived_inputs} function.
#'
#' @param config The config object that provides paths to required data
#' @return derived_input object
build_derived_inputs <- function(config){

  # NOTE: may fail
  tree <- load_tree_io(config@input@tree)

  species <- get_species_from_tree(tree)

  # NOTE: may fail
  queries <- read_queries(config@input@queries)

  # NOTE: may fail
  species_meta_list <- lapply(species, load_species, config@input)

  # NOTE: may fail
  synmap_meta_list <- lapply(species, load_synmaps, config@input@focal_species, config@input@synmaps)

  new(
    "derived_input",
    tree          = tree,
    focal_species = config@input@focal_species,
    queries       = queries,
    species       = species_meta_list,
    synmaps       = synmap_meta_list
  )

}
