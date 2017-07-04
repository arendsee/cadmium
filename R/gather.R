load_species <- function(species_name, input){
  # Primary data - required inputs to Fagin

  message(sprintf("Preprocessing data for %s ...", species_name))
  message(sprintf(" - loading genome", species_name))

  # TODO: handle errors
  dna_filename <- get_genome_filename ( species_name, input@fna_dir )
  dna <- load_dna(dna_filename)

  message(" - loading GFF")

  # TODO: handle errors
  # TODO: also, this is too slow, need to speed it up
  gff_filename <- get_gff_filename(species_name, input@gff_dir)
  gff <- load_gff(gff_filename)


  # Gene model independent data, derived only from genome
  message(" - finding N-strings")
  nstring <- derive_nstring(dna)

  message(" - finding finding genomic ORFs")
  orfgff <- derive_orfgff(dna)

  message(" - assembling proteins from genomic ORFs")
  orffaa <- derive_orffaa(dna, orfgff)

  # Gene model derived data
  message(" - assembling proteins from genome and GFF")
  aa <- derive_aa(dna, gff)

  message(" - assembling mRNA transcripts")
  trans <- derive_trans(dna, gff)

  message(" - finding ORFs on transcripts")
  transorfgff <- derive_transorfgff(trans)

  message(" - translating transcript ORFs")
  transorffaa <- derive_transorffaa(trans, transorfgff)

  message(" - summarizing data")

  message("   - gff.summary")
  gff.summary = summarize_gff(gff)

  message("   - dna.summary")
  dna.summary = summarize_dna(dna)

  message("   - aa.summary")
  aa.summary = summarize_faa(aa)

  message("   - trans.summary")
  trans.summary = summarize_dna(trans)

  message("   - orfgff.summary")
  orfgff.summary = summarize_granges(orfgff)

  message("   - orffaa.summary")
  orffaa.summary = summarize_faa(orffaa)

  message("   - transorfgff.summary")
  transorfgff.summary = summarize_granges(transorfgff)

  message("   - transorffaa.summary")
  transorffaa.summary = summarize_faa(transorffaa)

  message("   - nstring.summary")
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

  message(" - caching data")

  specfile <- new("species_data_files",
    gff.file         = to_cache( gff         , label="gff"         , group=species_name ),
    dna.file         = to_cache( dna         , label="dna"         , group=species_name ),
    aa.file          = to_cache( aa          , label="aa"          , group=species_name ),
    trans.file       = to_cache( trans       , label="trans"       , group=species_name ),
    orfgff.file      = to_cache( orfgff      , label="orfgff"      , group=species_name ),
    orffaa.file      = to_cache( orffaa      , label="orffaa"      , group=species_name ),
    transorfgff.file = to_cache( transorfgff , label="transorfgff" , group=species_name ),
    transorffaa.file = to_cache( transorffaa , label="transorffaa" , group=species_name ),
    nstring.file     = to_cache( nstring     , label="nstring"     , group=species_name ),
    specsum.file     = to_cache( specsum     , label="specsum"     , group=species_name )
  )

  new(
    "species_meta",
    files     = specfile,
    summaries = specsum
  )

}


load_species_with_cache <- function(species_name, ...){
  # Retrieve the cached data if present
  result <- from_cache(species_name, "species_meta")  
  if(is.null(result)){
    result <- load_species(species_name, ...) 
    to_cache(result, species_name, "species_meta")  
  }
  result
}


load_synmaps <- function(target_species, focal_species, syndir){

  # TODO: replace this hard-coded pattern
  filename <- get_synmap_filename(focal_species, target_species, syndir)

  synfile <- load_synmap(filename)

  synsum <- summarize_syn(synfile)

  new(
    "synteny_meta",
    synmap.file    = to_cache(synfile, label="synmap", group=target_species),
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
  species_meta_list <- BiocParallel::bplapply(species_names, load_species_with_cache, con@input)
  names(species_meta_list) <- species_names

  focal_species <- gsub(" ", "_", con@input@focal_species)

  # NOTE: may fail
  synmap_meta_list <- lapply(
    setdiff(species_names, focal_species),
    load_synmaps,
    focal_species = con@input@focal_species,
    syndir        = con@input@syn_dir
  )
  names(synmap_meta_list) <- setdiff(species_names, focal_species)

  new(
    "derived_input",
    tree          = tree,
    focal_species = focal_species,
    queries       = queries,
    species       = species_meta_list,
    synmaps       = synmap_meta_list
  )

}
