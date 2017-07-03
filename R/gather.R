load_species <- function(species_name, input){
  # Primary data - required inputs to Fagin

  cat(sprintf("Preprocessing data for %s ...\n", species_name))
  cat(sprintf(" - loading genome\n", species_name))

  # TODO: handle errors
  dna <- load_dna( get_genome_filename ( species_name, input@fna_dir ))

  cat(" - loading GFF\n")

  # TODO: handle errors
  # TODO: also, this is too slow, need to speed it up
  gff <- load_gff( get_gff_filename    ( species_name, input@gff_dir ))


  # Gene model independent data, derived only from genome
  cat(" - finding N-strings\n")
  nstring <- derive_nstring(dna)

  cat(" - finding finding genomic ORFs\n")
  orfgff <- derive_orfgff(dna)

  cat(" - assembling proteins from genomic ORFs\n")
  orffaa <- derive_orffaa(dna, orfgff)

  # Gene model derived data
  cat(" - assembling proteins from genome and GFF\n")
  aa <- derive_aa(dna, gff)

  cat(" - assembling mRNA transcripts\n")
  trans <- derive_trans(dna, gff)

  cat(" - finding ORFs on transcripts\n")
  transorfgff <- derive_transorfgff(trans)

  cat(" - translating transcript ORFs\n")
  transorffaa <- derive_transorffaa(trans, transorfgff)

  cat(" - summarizing data\n")

  cat("   - gff.summary\n")
  gff.summary = summarize_gff(gff)

  cat("   - dna.summary\n")
  dna.summary = summarize_dna(dna)

  cat("   - aa.summary\n")
  aa.summary = summarize_faa(aa)

  cat("   - trans.summary\n")
  trans.summary = summarize_dna(trans)

  cat("   - orfgff.summary\n")
  orfgff.summary = summarize_granges(orfgff)

  cat("   - orffaa.summary\n")
  orffaa.summary = summarize_faa(orffaa)

  cat("   - transorfgff.summary\n")
  transorfgff.summary = summarize_granges(transorfgff)

  cat("   - transorffaa.summary\n")
  transorffaa.summary = summarize_faa(transorffaa)

  cat("   - nstring.summary\n")
  nstring.summary = summarize_nstring(nstring)

  specsum <- new("species_summaries",
    gff.summary         = gff.summary,
    dna.summary         = dna.summary,
    aa.summary          = aa.summary,
    trans.summary       = trans.summary,
    orfgff.summary      = orfgff.summary,
    orffaa.summary      = orffaa.summary,
    transorfgff.summary = transorfgff.summary,
    transorffaa.summary = transorffaa.summary,
    nstring.summary     = nstring.summary
  )

  cat(" - caching data\n")

  specfile <- new("species_data_files",
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
  filename <- get_synmap_filename(focal_species, target_species, syndir)

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
#' @param con The config object that provides paths to required data
#' @return derived_input object
#' @export
load_data <- function(con){

  # NOTE: may fail
  tree <- load_tree(con@input@tree)

  # NOTE: may fail
  queries <- load_queries(con@input@query_gene_list)

  species_names <- tree$tip.label

  # NOTE: may fail
  species_meta_list <- lapply(species_names, load_species, con@input)

  # NOTE: may fail
  synmap_meta_list <- lapply(
    species_names,
    load_synmaps,
    target_species = con@input@focal_species,
    syndir         = con@input@synmaps
  )

  new(
    "derived_input",
    tree          = tree,
    focal_species = con@input@focal_species,
    queries       = queries,
    species       = species_meta_list,
    synmaps       = synmap_meta_list
  )

}
