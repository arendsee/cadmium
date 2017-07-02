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
    gff.file         = to_cache( gff         , label="gff"         , group=species_name ),
    dna.file         = to_cache( dna         , label="dna"         , group=species_name ),
    aa.file          = to_cache( aa          , label="aa"          , group=species_name ),
    trans.file       = to_cache( trans       , label="trans"       , group=species_name ),
    orfgff.file      = to_cache( orfgff      , label="orfgff"      , group=species_name ),
    orffaa.file      = to_cache( orffaa      , label="orffaa"      , group=species_name ),
    transorfgff.file = to_cache( transorfgff , label="transorfgff" , group=species_name ),
    transorffaa.file = to_cache( transorffaa , label="transorffaa" , group=species_name ),
    nstring.file     = to_cache( nstring     , label="nstring"     , group=species_name )
  )

  new(
    "species_meta",
    files     = specfile,
    summaries = specsum
  )

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
  tree <- load_tree(config@input@tree)

  # NOTE: may fail
  queries <- load_queries(config@input@queries)

  species_names <- tree$tip.label

  # NOTE: may fail
  species_meta_list <- lapply(species_names, load_species, config@input)

  # NOTE: may fail
  synmap_meta_list <- lapply(
    species_names,
    load_synmaps,
    taget_species = config@input@focal_species,
    syndir        = config@input@synmaps
  )

  new(
    "derived_input",
    tree          = tree,
    focal_species = config@input@focal_species,
    queries       = queries,
    species       = species_meta_list,
    synmaps       = synmap_meta_list
  )

}
